classdef MMn < MomentModel
    %MMn: Mixed-moment minimum-entropy models of arbitrary order
    
    
    properties
        ChebyshevQuadMu
        SphericalHarmonicsQuadMu
        alpha
        
    end
    
    methods
        
        function obj = MMn(varargin)
            
            %Inputparser
            p = inputParser;
            p.KeepUnmatched=true;
            p.addParameter('QuadratureOrder',20);
            p.parse(varargin{:});
            
            Fields = p.Unmatched;
            
            order = p.Results.QuadratureOrder;
            nqs = (order+1+double(mod(order,2)==0))/2+1;
            load HalfRangeChebyshevPolys
            nqs = min(nqs,length(xwd)); %#ok<USENS>
            x = xwd{nqs}(:,1);
            w = xwd{nqs}(:,2);
            
            [Qpp,Qmp,Qmm,Qpm] = SphericalQuadratureQM(order);
            Q = [Qpp,Qmp,Qmm,Qpm];
            Fields.QuadraturePoints = Q(1:2,:); %Overwrite given quadrature rule
            Fields.QuadratureWeights = Q(3,:);
            
            obj = obj@MomentModel(Fields);
            obj.QuadratureOrder = order;
            obj.positivityFlag = 1;
            obj.ChebyshevQuadMu.mu = x;
            obj.ChebyshevQuadMu.wmu = 2*w;
            obj.SphericalHarmonicsQuadMu.mu = sqrt(1-x.^2);
            obj.SphericalHarmonicsQuadMu.wmu = 2*w.*x;
            
        end
        
        
        function [ BasisFunctions,BasisFunctionsAtQuadrature,BasisFunctionsAtQuadratureSqueezed] = generateBasisFunctions(obj)
            
            BasisFunctions = [];
            %
            mu = permute(obj.QuadraturePoints(1,:),[1 3 2]);
            phi = permute(obj.QuadraturePoints(2,:),[1 3 2]);
            if min(phi(:))<0
                phi = phi+pi;
            end
            lmax = obj.MomentOrder;
            
            Omegax = sqrt(1-mu.^2).*cos(phi);
            Omegay = sqrt(1-mu.^2).*sin(phi);
            
            indpp = phi<=pi/2;
            indmp = phi>=pi/2 & phi<=pi;
            indmm = phi>=pi & phi<=3*pi/2;
            indpm = phi>=3*pi/2 & phi<=2*pi;
            indSxp = indpp | indpm;
            indSxm = indmp | indmm;
            indSyp = indpp | indmp;
            indSym = indpm | indmm;
            
            BFSxp = bsxfun(@power,Omegax.*indSxp,transpose(1:lmax));
            BFSxm = bsxfun(@power,Omegax.*indSxm,transpose(1:lmax));
            BFSyp = bsxfun(@power,Omegay.*indSyp,transpose(1:lmax));
            BFSym = bsxfun(@power,Omegay.*indSym,transpose(1:lmax));
            [I,J] = meshgrid(1:lmax-1);
            ind = I+J<=lmax;
            I = I(ind);
            J = J(ind);
            if lmax>1
                BFpp = bsxfun(@power,Omegax.*indpp,I).*bsxfun(@power,Omegay.*indpp,J);
                BFmp = bsxfun(@power,Omegax.*indmp,I).*bsxfun(@power,Omegay.*indmp,J);
                BFmm = bsxfun(@power,Omegax.*indmm,I).*bsxfun(@power,Omegay.*indmm,J);
                BFpm = bsxfun(@power,Omegax.*indpm,I).*bsxfun(@power,Omegay.*indpm,J);
                
                BasisFunctionsAtQuadrature = cat(1,1+0*Omegax,BFSxp,BFSxm,BFSyp,BFSym,BFpp,BFmp,BFmm,BFpm);
            else
                BasisFunctionsAtQuadrature = cat(1,1+0*Omegax,BFSxp,BFSxm,BFSyp,BFSym);
            end
            BasisFunctionsAtQuadratureSqueezed = squeeze(BasisFunctionsAtQuadrature);
        end
        
        
        function [psi,out] = ReconstructDistribution(obj,u)
            persistent beta BM
            
            beta_size = size(u);
            if isempty(beta) || ~isequal(size(beta),beta_size)
                beta = repmat([log(1/4/pi);zeros(size(u,1)-1,1)],[1,beta_size(2)]);
                BM = repmat(eye(size(u,1)),[1 1 beta_size(2)]);
            end
            
            
            OptimPara = obj.OptimPara;
            psi = zeros([size(obj.BasisFunctionsAtQuadrature,3),size(u,2)]);
            
            for i=1:size(u,2);
                U = u(:,i);
                [alphaF,beta(:,i),BM(:,:,i)] = dualAdaptivePoly(U, beta(:,i),BM(:,:,i),OptimPara);
                psi(:,i) = U(1)*exp(alphaF'*obj.BasisFunctionsAtQuadratureSqueezed);
            end
            
            out = [];
        end
        
        function [UC,J] = NeutronKernel(obj,~,U)
            UC = bsxfun(@times,obj.phiIso/obj.phiIso(1),U(1,:,:,:))-U;
            J = repmat([obj.phiIso/obj.phiIso(1),zeros(obj.NumberMoments,obj.NumberMoments-1)]+diag(-ones(obj.NumberMoments,1)),[1 1 size(UC,2)]);
        end
        
        
        
        function [LB,J] = LaplaceBeltrami2(obj,~,U) %Linearization
            persistent S
            
            if isempty(S) || ~isequal(size(S,1),size(U,1))
                [~, ~,S] = generateMPnClosure(obj.MomentOrder);
            end
            LB = permute(sum(bsxfun(@times,S,permute(U,[5 1 2 3 4])),2),[1 3 4 5 2]);
            J = repmat(S,[1 1 size(U,2)]);
        end
        %
        function [LB,Jc] = LaplaceBeltrami(obj,psi,U)
            persistent DT DT2
                
            if obj.MomentOrder > 3
                [LB,Jc] = LaplaceBeltrami2(obj,psi,U); %Solving the nonlinear system is too expensive -> MPN
            else
                
                %New quadrature rules and interpolation
                Jc = []; %Let the nonlinear solver calculate the derivatives for us
                mu = obj.QuadraturePoints(1,:);
                phi = obj.QuadraturePoints(2,:);
                
                mu2 = obj.ChebyshevQuadMu.mu;
                w = obj.ChebyshevQuadMu.wmu;
                I = mu==min(mu);
                mu = [mu,0*mu(I)];
                phi = [phi,phi(I)];
                psival = [psi;psi(I,:)];
                xyz = [sqrt(1-mu.^2).*cos(phi);sqrt(1-mu.^2).*sin(phi);mu];
                
                if isempty(DT)
                    DT  = delaunayTriangulation(xyz(1:2,:)');
                end
                phi2 = 0*mu2;
                xyz_0 = [sqrt(1-mu2.^2).*cos(phi2),sqrt(1-mu2.^2).*sin(phi2),mu2];
                [ti,bc] = pointLocation(DT,xyz_0(:,1:2));
                trival = reshape(psival(DT.ConnectivityList(ti,:),:),[size(bc),size(psival,2)]);
                psival1 = permute(sum(bsxfun(@times,bc,trival),2),[1 3 2]); %psi(mu,0)
                
                phi2 = pi/2+0*mu2;
                xyz_0 = [sqrt(1-mu2.^2).*cos(phi2),sqrt(1-mu2.^2).*sin(phi2),mu2];
                [ti,bc] = pointLocation(DT,xyz_0(:,1:2));
                trival = reshape(psival(DT.ConnectivityList(ti,:),:),[size(bc),size(psival,2)]);
                psival2 = permute(sum(bsxfun(@times,bc,trival),2),[1 3 2]); %psi(mu,pi/2)
                
                phi2 = pi+0*mu2;
                xyz_0 = [sqrt(1-mu2.^2).*cos(phi2),sqrt(1-mu2.^2).*sin(phi2),mu2];
                [ti,bc] = pointLocation(DT,xyz_0(:,1:2));
                trival = reshape(psival(DT.ConnectivityList(ti,:),:),[size(bc),size(psival,2)]);
                psival3 = permute(sum(bsxfun(@times,bc,trival),2),[1 3 2]); %psi(mu,pi)
                
                phi2 = 3*pi/2+0*mu2;
                xyz_0 = [sqrt(1-mu2.^2).*cos(phi2),sqrt(1-mu2.^2).*sin(phi2),mu2];
                [ti,bc] = pointLocation(DT,xyz_0(:,1:2));
                trival = reshape(psival(DT.ConnectivityList(ti,:),:),[size(bc),size(psival,2)]);
                psival4 = permute(sum(bsxfun(@times,bc,trival),2),[1 3 2]); %psi(mu,3pi/2)
                
                %with chebyshev quadrature rule (weight function is
                %1/sqrt(1-mu^2)
                int1 = sum(bsxfun(@times,w,psival1),1);
                int2 = sum(bsxfun(@times,w,psival2),1);
                int3 = sum(bsxfun(@times,w,psival3),1);
                int4 = sum(bsxfun(@times,w,psival4),1);
                
                mu = obj.QuadraturePoints(1,:);
                phi =obj.QuadraturePoints(2,:);
                if min(phi(:))<0
                    phi = phi+pi;
                end
                Omegax = sqrt(1-mu.^2).*cos(phi);
                Omegay = sqrt(1-mu.^2).*sin(phi);
                indpp = phi<=pi/2;
                indmp = phi>pi/2 & phi<=pi;
                indmm = phi>pi & phi<=3*pi/2;
                indpm = phi>3*pi/2 & phi<=2*pi;
                
                Basis = [double(indpp);double(indmp);double(indmm);double(indpm)];
                BasisHalf = [double(indpp | indpm);double(indmp | indmm);double(indpp | indmp);double(indpm | indmm)];
                
                LB = 0*U; %initializes u0
                lmax = obj.MomentOrder;
                l=1;
                LB(l+1,:) = -l*(l+1)*U(l+1,:) + int2+int4; %Sxp
                LB(lmax+l+1,:) = -l*(l+1)*U(lmax+l+1,:)-int2-int4; %Sxm
                LB(2*lmax+l+1,:) = -l*(l+1)*U(2*lmax+l+1,:)+int1+int3; %Syp
                LB(3*lmax+l+1,:) = -l*(l+1)*U(3*lmax+l+1,:)-int1-int3; %Sym
                
                if lmax>1
                    l=2;
                    u = projectBasis(psi,obj.QuadratureWeightsSqueezed,BasisHalf,3);
                    LB(l+1,:) = -l*(l+1)*U(l+1,:)+l*(l-1)*u(1,:); %Sxp
                    LB(lmax+l+1,:) = -l*(l+1)*U(lmax+l+1,:)+l*(l-1)*u(2,:); %Sxm
                    LB(2*lmax+l+1,:) = -l*(l+1)*U(2*lmax+l+1,:)+l*(l-1)*u(3,:); %Syp
                    LB(3*lmax+l+1,:) = -l*(l+1)*U(3*lmax+l+1,:)+l*(l-1)*u(4,:); %Sym
                end
                
                for l=3:lmax
                    LB(l+1,:) = -l*(l+1)*U(l+1,:)+l*(l-1)*U(l-1,:); %Sxp
                    LB(lmax+l+1,:) = -l*(l+1)*U(lmax+l+1,:)+l*(l-1)*U(lmax+l-1,:); %Sxm
                    LB(2*lmax+l+1,:) = -l*(l+1)*U(2*lmax+l+1,:)+l*(l-1)*U(2*lmax+l-1,:); %Syp
                    LB(3*lmax+l+1,:) = -l*(l+1)*U(3*lmax+l+1,:)+l*(l-1)*U(3*lmax+l-1,:); %Sym
                end
                
                
                mu2 = obj.SphericalHarmonicsQuadMu.mu; %Switch to standard quadrature since Chebyshev has no advantage for weight function 1/(1-mu^2)
                w = obj.SphericalHarmonicsQuadMu.wmu; 
                
                I = mu==min(mu);
                mu = [mu,0*mu(I)];
                phi = [phi,phi(I)];
                psival = [psi;psi(I,:)];
                xyz = [sqrt(1-mu.^2).*cos(phi);sqrt(1-mu.^2).*sin(phi);mu];
                
                if isempty(DT2)
                    DT2  = delaunayTriangulation(xyz(1:2,:)');
                end
                phi2 = 0*mu2;
                xyz_0 = [sqrt(1-mu2.^2).*cos(phi2),sqrt(1-mu2.^2).*sin(phi2),mu2];
                [ti,bc] = pointLocation(DT2,xyz_0(:,1:2));
                trival = reshape(psival(DT2.ConnectivityList(ti,:),:),[size(bc),size(psival,2)]);
                psival1 = permute(sum(bsxfun(@times,bc,trival),2),[1 3 2]); %psi(mu,0)
                
                phi2 = pi/2+0*mu2;
                xyz_0 = [sqrt(1-mu2.^2).*cos(phi2),sqrt(1-mu2.^2).*sin(phi2),mu2];
                [ti,bc] = pointLocation(DT2,xyz_0(:,1:2));
                trival = reshape(psival(DT2.ConnectivityList(ti,:),:),[size(bc),size(psival,2)]);
                psival2 = permute(sum(bsxfun(@times,bc,trival),2),[1 3 2]); %psi(mu,pi/2)
                
                phi2 = pi+0*mu2;
                xyz_0 = [sqrt(1-mu2.^2).*cos(phi2),sqrt(1-mu2.^2).*sin(phi2),mu2];
                [ti,bc] = pointLocation(DT2,xyz_0(:,1:2));
                trival = reshape(psival(DT2.ConnectivityList(ti,:),:),[size(bc),size(psival,2)]);
                psival3 = permute(sum(bsxfun(@times,bc,trival),2),[1 3 2]); %psi(mu,pi)
                
                phi2 = 3*pi/2+0*mu2;
                xyz_0 = [sqrt(1-mu2.^2).*cos(phi2),sqrt(1-mu2.^2).*sin(phi2),mu2];
                [ti,bc] = pointLocation(DT2,xyz_0(:,1:2));
                trival = reshape(psival(DT2.ConnectivityList(ti,:),:),[size(bc),size(psival,2)]);
                psival4 = permute(sum(bsxfun(@times,bc,trival),2),[1 3 2]); %psi(mu,3pi/2)
                
                
                Omegax1 = sqrt(1-mu2.^2)'.*cos(0);
                Omegay1 = sqrt(1-mu2.^2)'.*sin(0);
                Omegax2 = sqrt(1-mu2.^2)'.*cos(pi/2);
                Omegay2 = sqrt(1-mu2.^2)'.*sin(pi/2);
                Omegax3 = sqrt(1-mu2.^2)'.*cos(pi);
                Omegay3 = sqrt(1-mu2.^2)'.*sin(pi);
                Omegax4 = sqrt(1-mu2.^2)'.*cos(3*pi/2);
                Omegay4 = sqrt(1-mu2.^2)'.*sin(3*pi/2);
                
                [I,J] = meshgrid(1:lmax-1); %Quarter moments
                ind = I+J<=lmax;
                I = I(ind);
                J = J(ind);
                s = size(I,1);
                
                for i=1:s
                    m = I(i)+J(i);
                    l = I(i);
                    
                    %Standard quadrature
                    int1 = sum(bsxfun(@times,w.*transpose(Omegax1.^(l+1).*Omegay1.^(m-l-1)*(l-m)+l*Omegax1.^(l-1).*Omegay1.^(m-l+1))./(1-mu2.^2),psival1),1); 
                    int2 = sum(bsxfun(@times,w.*transpose(Omegax2.^(l+1).*Omegay2.^(m-l-1)*(l-m)+l*Omegax2.^(l-1).*Omegay2.^(m-l+1))./(1-mu2.^2),psival2),1);
                    int3 = sum(bsxfun(@times,w.*transpose(Omegax3.^(l+1).*Omegay3.^(m-l-1)*(l-m)+l*Omegax3.^(l-1).*Omegay3.^(m-l+1))./(1-mu2.^2),psival3),1);
                    int4 = sum(bsxfun(@times,w.*transpose(Omegax4.^(l+1).*Omegay4.^(m-l-1)*(l-m)+l*Omegax4.^(l-1).*Omegay4.^(m-l+1))./(1-mu2.^2),psival4),1);
                    
                    LB(1+4*lmax+i,:) = -U(1+4*lmax+i,:)*(m+m^2)+int2-int1;
                    LB(1+4*lmax+i+s,:) = -U(1+4*lmax+i+s,:)*(m+m^2)+int3-int2;
                    LB(1+4*lmax+i+2*s,:) = -U(1+4*lmax+i+2*s,:)*(m+m^2)+int4-int3;
                    LB(1+4*lmax+i+3*s,:) = -U(1+4*lmax+i+3*s,:)*(m+m^2)+int1-int4;
                    if m-l>1
                        u = projectBasis(bsxfun(@times,psi,transpose(Omegax.^l.*Omegay.^(m-l-2))),obj.QuadratureWeightsSqueezed,Basis,3);
                        LB(1+4*lmax+i,:) = LB(1+4*lmax+i,:)+(m-l)*(m-l-1)*u(1,:);
                        LB(1+4*lmax+i+s,:) = LB(1+4*lmax+i+s,:)+(m-l)*(m-l-1)*u(2,:);
                        LB(1+4*lmax+i+2*s,:) = LB(1+4*lmax+i+2*s,:)+(m-l)*(m-l-1)*u(3,:);
                        LB(1+4*lmax+i+3*s,:) = LB(1+4*lmax+i+3*s,:)+(m-l)*(m-l-1)*u(4,:);
                    end
                    if l>1
                        u = projectBasis(bsxfun(@times,psi,transpose(Omegax.^(l-2).*Omegay.^(m-l))),obj.QuadratureWeightsSqueezed,Basis,3);
                        LB(1+4*lmax+i,:) = LB(1+4*lmax+i,:)+l*(l-1)*u(1,:);
                        LB(1+4*lmax+i+s,:) = LB(1+4*lmax+i+s,:)+l*(l-1)*u(2,:);
                        LB(1+4*lmax+i+2*s,:) = LB(1+4*lmax+i+2*s,:)+l*(l-1)*u(3,:);
                        LB(1+4*lmax+i+3*s,:) = LB(1+4*lmax+i+3*s,:)+l*(l-1)*u(4,:);
                    end
                end
                
            end                       
        end                               
    end
    
end

