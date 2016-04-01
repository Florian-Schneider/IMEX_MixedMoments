classdef MK1 < MomentModel
    %MK1: Mixed-moment Kershaw closure of first order with linearized
    %Laplace-Beltrami operator
    
    properties
    end
    
    methods
        
        function obj = MK1(varargin)
            
            %Inputparser
            p = inputParser;
            p.KeepUnmatched=true;
            p.addParameter('QuadratureOrder',20);
            p.parse(varargin{:});
            
            Fields = p.Unmatched;
            
            order = p.Results.QuadratureOrder;
                        
            [Qpp,Qmp,Qmm,Qpm] = SphericalQuadratureQM(order);
            Q = [Qpp,Qmp,Qmm,Qpm];
            Fields.QuadraturePoints = Q(1:2,:); %Overwrite given quadrature rule
            Fields.QuadratureWeights = Q(3,:);
            Fields.MomentOrder = 1; %Force moment order.
            
            obj = obj@MomentModel(Fields);
            obj.positivityFlag = 1;
            
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
        
        function [FXP,FXM,FYP,FYM,Uc,Jc] = Closure(obj,U,R,scatteringtrue)
            tol = 1e-10;
            
            xi = -(2*pi-3*pi*sqrt((2))+4)/(3*pi*(sqrt((2))-1));
            a = 1;
            b = (-(2)^(1/2)*(a - xi));
            
            
            E = U(1,:,:,:);
            Fxp = U(2,:,:,:);
            Fxm = U(3,:,:,:);
            Fyp = U(4,:,:,:);
            Fym = U(5,:,:,:);
            
            fxp = Fxp./E;
            fxm = Fxm./E;
            fyp = Fyp./E;
            fym = Fym./E;
            
            fppx = (fxp-fxm);
            fppy = (fyp-fym);
            fpmx = (fxp-fxm);
            fpmy = -(fyp-fym);
            fmpx = -(fxp-fxm);
            fmpy = (fyp-fym);
            fmmx = -(fxp-fxm);
            fmmy = -(fyp-fym);
            
            a1pp = (fppx.*fppy)*(16/(3*pi) - 8/3).*(fppx.^2+fppy.^2 - 1);
            n = sqrt(fppx.^2+fppy.^2);
            zeta = a+b*n;
            a2pp = -a1pp + zeta.*n.^2+(1-zeta).*n;
            
            nxx = fppx.*fppx./(fppx.^2+fppy.^2);
            nxy = fppx.*fppy./(fppx.^2+fppy.^2);
            nyy = fppy.*fppy./(fppx.^2+fppy.^2);
            nxx((fppx.^2+fppy.^2).^2<=tol) = 0;
            nxy((fppx.^2+fppy.^2).^2<=tol) = 0;
            nyy((fppx.^2+fppy.^2).^2<=tol) = 0;
            Dppxx = a1pp+a2pp.*nxx;
            Dppxy = a2pp.*nxy;
            Dppyy = a1pp+a2pp.*nyy;
            
            a1mp = -(fmpx.*fmpy)*(16/(3*pi) - 8/3).*(fmpx.^2+fmpy.^2 - 1);
            n = sqrt(fmpx.^2+fmpy.^2);
            zeta = a+b*n;
            a2mp = -a1mp + zeta.*n.^2+(1-zeta).*n;
            
            nxx = fmpx.*fmpx./(fmpx.^2+fmpy.^2);
            nxy = fmpx.*fmpy./(fmpx.^2+fmpy.^2);
            nyy = fmpy.*fmpy./(fmpx.^2+fmpy.^2);
            nxx((fmpx.^2+fmpy.^2).^2<=tol) = 0;
            nxy((fmpx.^2+fmpy.^2).^2<=tol) = 0;
            nyy((fmpx.^2+fmpy.^2).^2<=tol) = 0;
            Dmpxx = a1mp+a2mp.*nxx;
            Dmpxy = a2mp.*nxy;
            Dmpyy = a1mp+a2mp.*nyy;
            
            a1mm = (fmmx.*fmmy)*(16/(3*pi) - 8/3).*(fmmx.^2+fmmy.^2 - 1);
            n = sqrt(fmmx.^2+fmmy.^2);
            zeta = a+b*n;
            a2mm = -a1mm + zeta.*n.^2+(1-zeta).*n;
            
            nxx = fmmx.*fmmx./(fmmx.^2+fmmy.^2);
            nxy = fmmx.*fmmy./(fmmx.^2+fmmy.^2);
            nyy = fmmy.*fmmy./(fmmx.^2+fmmy.^2);
            nxx((fmmx.^2+fmmy.^2).^2<=tol) = 0;
            nxy((fmmx.^2+fmmy.^2).^2<=tol) = 0;
            nyy((fmmx.^2+fmmy.^2).^2<=tol) = 0;
            Dmmxx = a1mm+a2mm.*nxx;
            Dmmxy = a2mm.*nxy;
            Dmmyy = a1mm+a2mm.*nyy;
            
            a1pm = -(fpmx.*fpmy)*(16/(3*pi) - 8/3).*(fpmx.^2+fpmy.^2 - 1);
            n = sqrt(fpmx.^2+fpmy.^2);
            zeta = a+b*n;
            a2pm = -a1pm + zeta.*n.^2+(1-zeta).*n;
            
            nxx = fpmx.*fpmx./(fpmx.^2+fpmy.^2);
            nxy = fpmx.*fpmy./(fpmx.^2+fpmy.^2);
            nyy = fpmy.*fpmy./(fpmx.^2+fpmy.^2);
            nxx((fpmx.^2+fpmy.^2).^2<=tol) = 0;
            nxy((fpmx.^2+fpmy.^2).^2<=tol) = 0;
            nyy((fpmx.^2+fpmy.^2).^2<=tol) = 0;
            Dpmxx = a1pm+a2pm.*nxx;
            Dpmxy = a2pm.*nxy;
            Dpmyy = a1pm+a2pm.*nyy;
            
            
            lambda1 = fyp./(fyp-fym);
            lambda1(fyp-fym<tol) = 0;
            lambda1 = max(min(lambda1,1),0);
            lambda2 = fxp./(fxp-fxm);
            lambda2(fxp-fxm<tol) = 0;
            lambda2 = max(min(lambda2,1),0);
            
            Pxp = E.*(lambda1.*lambda2.*Dppxx + (1-lambda1).*lambda2.*Dpmxx);
            Pxm = E.*(lambda1.*(1-lambda2).*Dmpxx + (1-lambda1).*(1-lambda2).*Dmmxx);
            Pyp = E.*(lambda1.*lambda2.*Dppyy + lambda1.*(1-lambda2).*Dmpyy);
            Pym = E.*((1-lambda1).*lambda2.*Dpmyy + (1-lambda1).*(1-lambda2).*Dmmyy);
            Ppp = E.*lambda2.*lambda1.*Dppxy;
            Pmp = E.*(1-lambda2).*lambda1.*Dmpxy;
            Pmm = E.*(1-lambda2).*(1-lambda1).*Dmmxy;
            Ppm = E.*lambda2.*(1-lambda1).*Dpmxy;
            
            FXP = cat(1,Fxp+0*Fxm,Pxp,0*Pxm,Ppp+0*Pmp,Ppm+0*Pmm);
            FXM = cat(1,0*Fxp+Fxm,0*Pxp,Pxm,0*Ppp+Pmp,0*Ppm+Pmm);
            FYP = cat(1,Fyp+0*Fym,Ppp+0*Ppm,Pmp+0*Pmm,Pyp,0*Pym);
            FYM = cat(1,0*Fyp+Fym,0*Ppp+Ppm,0*Pmp+Pmm,0*Pyp,Pym);
            if scatteringtrue
                [Uc,Jc] = R.Collision([],U);
                Uc = reshape(Uc,size(U));
            else
                Uc = 0*U;
                Jc = zeros(obj.NumberMoments,obj.NumberMoments,size(U,2));
            end
            
        end
        
        
        
        function [psi] = ReconstructDistribution(obj,u)
            error('does not exist for Kershaw models')
            
        end
        
    
        
        function [UC,J] = NeutronKernel(obj,~,U)
            UC = bsxfun(@times,obj.phiIso/obj.phiIso(1),U(1,:,:,:))-U;
            J = repmat([obj.phiIso/obj.phiIso(1),zeros(obj.NumberMoments,obj.NumberMoments-1)]+diag(-ones(obj.NumberMoments,1)),[1 1 size(UC,2)]);
        end
        
        
        
        function [LB,J] = LaplaceBeltrami(obj,~,U)
            persistent S
            
            if isempty(S) || ~isequal(size(S,1),size(U,1))
                [~, ~,S] = generateMPnClosure(obj.MomentOrder);
            end
            LB = permute(sum(bsxfun(@times,S,permute(U,[5 1 2 3 4])),2),[1 3 4 5 2]);
            J = repmat(S,[1 1 size(U,2)]); %Jacobian
        end
        
       
        
        
    end
    
end

