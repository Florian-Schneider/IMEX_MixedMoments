classdef QMn
    
    
    properties
        BasisFunctions
        BasisFunctionsAtQuadrature
        BasisFunctionsAtQuadratureSqueezed
        MomentOrder
        NumberMoments
        nq
        OptimParapp
        OptimParamp
        OptimParamm
        OptimParapm
        Omegaxpp
        Omegaypp
        Omegaxmp
        Omegaymp
        Omegaxmm
        Omegaymm
        Omegaxpm
        Omegaypm
        phiIso
        positivityFlag
        Density
        
        QuadraturePoints_pp
        QuadratureWeights_pp
        QuadraturePoints_mp
        QuadratureWeights_mp
        QuadraturePoints_pm
        QuadratureWeights_pm
        QuadraturePoints_mm
        QuadratureWeights_mm
        QuadraturePoints
        QuadratureWeights
        
        TimeStepping
    end
    
    methods
        
        function obj = QMn(varargin)
            
            %Inputparser
            p = inputParser;
            p.addParameter('MomentOrder',[]);
            p.addParameter('QuadratureOrder',20);
            p.KeepUnmatched=true;
            p.parse(varargin{:});
            
            
            [Qpp,Qmp,Qmm,Qpm] = SphericalQuadratureQM(p.Results.QuadratureOrder);
            
            obj.positivityFlag = 1;
            
            NN = obj.NumberMoments/4;
            
            obj.QuadraturePoints_pp = Qpp(1:2,:);
            obj.QuadratureWeights_pp = Qpp(3,:);
            obj.QuadraturePoints_mp = Qmp(1:2,:);
            obj.QuadratureWeights_mp = Qmp(3,:);
            obj.QuadraturePoints_mm = Qmm(1:2,:);
            obj.QuadratureWeights_mm = Qmm(3,:);
            obj.QuadraturePoints_pm = Qpm(1:2,:);
            obj.QuadratureWeights_pm = Qpm(3,:);
            obj.MomentOrder=p.Results.MomentOrder;
            Q = [Qpp,Qmp,Qmm,Qpm];
            obj.QuadraturePoints = Q(1:2,:); %Overwrite given quadrature rule
            obj.QuadratureWeights = Q(3,:);
            
            
            obj.nq = 4*size(Qpp,2);
            [ obj.BasisFunctions,obj.BasisFunctionsAtQuadrature,obj.BasisFunctionsAtQuadratureSqueezed] = obj.generateBasisFunctions;
            
            obj.NumberMoments = size(obj.BasisFunctionsAtQuadrature,1);
            
            obj.Omegaxpp = transpose(sqrt(1-obj.QuadraturePoints_pp(1,:).^2).*cos(obj.QuadraturePoints_pp(2,:)));
            obj.Omegaypp = transpose(sqrt(1-obj.QuadraturePoints_pp(1,:).^2).*sin(obj.QuadraturePoints_pp(2,:)));
            obj.Omegaxmp = transpose(sqrt(1-obj.QuadraturePoints_mp(1,:).^2).*cos(obj.QuadraturePoints_mp(2,:)));
            obj.Omegaymp = transpose(sqrt(1-obj.QuadraturePoints_mp(1,:).^2).*sin(obj.QuadraturePoints_mp(2,:)));
            obj.Omegaxmm = transpose(sqrt(1-obj.QuadraturePoints_mm(1,:).^2).*cos(obj.QuadraturePoints_mm(2,:)));
            obj.Omegaymm = transpose(sqrt(1-obj.QuadraturePoints_mm(1,:).^2).*sin(obj.QuadraturePoints_mm(2,:)));
            obj.Omegaxpm = transpose(sqrt(1-obj.QuadraturePoints_pm(1,:).^2).*cos(obj.QuadraturePoints_pm(2,:)));
            obj.Omegaypm = transpose(sqrt(1-obj.QuadraturePoints_pm(1,:).^2).*sin(obj.QuadraturePoints_pm(2,:)));
            
            
            psi = ones(obj.nq,1);
            obj.phiIso = obj.ProjectToBasis(psi);
            [obj.OptimParapp,obj.OptimParamp,obj.OptimParamm,obj.OptimParapm] = generateOptimPara(obj);
            obj.positivityFlag = 1;
            sizecap = @(s) s(2:end);
            obj.Density = @(U) reshape(U(1,:)+U(obj.NumberMoments/4+1,:)+U(2*obj.NumberMoments/4+1,:)+U(3*obj.NumberMoments/4+1,:),[1,sizecap(size(U))]);
        end
        
        
        function [ BasisFunctions,BasisFunctionsAtQuadrature,BasisFunctionsAtQuadratureSqueezed] = generateBasisFunctions(obj)
            
            BasisFunctions = [];
            
            mupp = permute(obj.QuadraturePoints_pp(1,:),[1 3 2]);
            phipp = permute(obj.QuadraturePoints_pp(2,:),[1 3 2]);
            mump = permute(obj.QuadraturePoints_mp(1,:),[1 3 2]);
            phimp = permute(obj.QuadraturePoints_mp(2,:),[1 3 2]);
            mumm = permute(obj.QuadraturePoints_mm(1,:),[1 3 2]);
            phimm = permute(obj.QuadraturePoints_mm(2,:),[1 3 2]);
            mupm = permute(obj.QuadraturePoints_pm(1,:),[1 3 2]);
            phipm = permute(obj.QuadraturePoints_pm(2,:),[1 3 2]);
            
            lmax = obj.MomentOrder;
            
            Omegaxpp = sqrt(1-mupp.^2).*cos(phipp); %#ok<*PROP>
            Omegaypp = sqrt(1-mupp.^2).*sin(phipp);
            Omegaxmp = sqrt(1-mump.^2).*cos(phimp);
            Omegaymp = sqrt(1-mump.^2).*sin(phimp);
            Omegaxmm = sqrt(1-mumm.^2).*cos(phimm);
            Omegaymm = sqrt(1-mumm.^2).*sin(phimm);
            Omegaxpm = sqrt(1-mupm.^2).*cos(phipm);
            Omegaypm = sqrt(1-mupm.^2).*sin(phipm);
            
            cnt = 1;
            
            BFpp = zeros(2*lmax+1,1,length(Omegaxpp));
            BFmp = BFpp;
            BFmm = BFpp;
            BFpm = BFpp;
            for l=0:lmax
                for m=l:-1:0
                    BFpp(cnt,1,:) = Omegaxpp.^m.*Omegaypp.^(l-m);
                    BFmp(cnt,1,:) = Omegaxmp.^m.*Omegaymp.^(l-m);
                    BFmm(cnt,1,:) = Omegaxmm.^m.*Omegaymm.^(l-m);
                    BFpm(cnt,1,:) = Omegaxpm.^m.*Omegaypm.^(l-m);
                    cnt = cnt+1;
                end
            end
            
            
            
            BasisFunctionsAtQuadrature = cat(1,BFpp,BFmp,BFmm,BFpm);
            BasisFunctionsAtQuadratureSqueezed = squeeze(BasisFunctionsAtQuadrature);
        end
        
        
        function u0 = ProjectToDensity(obj,U,c)
            %Moments are stored in the c-th component of U. Function
            %projects to u0
            u0 = obj.Density(U);
        end
        
        function [OptimParapp,OptimParamp,OptimParamm,OptimParapm] = generateOptimPara(obj)
            NN = obj.NumberMoments/4;
            OptimPara.itermax = 300;
            OptimPara.rList = [0,1e-4,1e-2,1e-1];
            OptimPara.tol = 1e-9; % tolerance on the norm of the gradient in the
            %                       stopping criterion
            OptimPara.rList = [0, 1e-8, 1e-6, 1e-4];% series of regularization values to try
            OptimPara.k0 = 50; % number of Newton iterations to try before increasing r
            OptimPara.itermax = max(300, OptimPara.k0 * length(OptimPara.rList)); % maximum number of iterations
            OptimPara.rMax = max(OptimPara.rList);
            r0 = 1e-4; % should be positive and fairly small
            kappa = 0.5; % should be in ]0, 1[
            ellMax = ceil(log(eps / r0) / log(kappa));
            OptimPara.rList2 = [r0 * kappa.^(0:ellMax), 0];
            OptimParapp = OptimPara;
            OptimParamp = OptimPara;
            OptimParamm = OptimPara;
            OptimParapm = OptimPara;
            
            OptimParapp.p = obj.BasisFunctionsAtQuadratureSqueezed(1:NN,:);
            OptimParamp.p = obj.BasisFunctionsAtQuadratureSqueezed(NN+1:2*NN,:);
            OptimParamm.p = obj.BasisFunctionsAtQuadratureSqueezed(2*NN+1:3*NN,:);
            OptimParapm.p = obj.BasisFunctionsAtQuadratureSqueezed(3*NN+1:4*NN,:);
            
            OptimParapp.phiIso = obj.phiIso(1:NN,:);
            OptimParamp.phiIso = obj.phiIso(NN+1:2*NN,:);
            OptimParamm.phiIso = obj.phiIso(2*NN+1:3*NN,:);
            OptimParapm.phiIso = obj.phiIso(3*NN+1:4*NN,:);
            
            
            OptimParapp.w = obj.QuadratureWeights_pp;
            OptimParamp.w = obj.QuadratureWeights_mp;
            OptimParamm.w = obj.QuadratureWeights_mm;
            OptimParapm.w = obj.QuadratureWeights_pm;
            
            OptimParapp.NN = NN;
            OptimParamp.NN = NN;
            OptimParamm.NN = NN;
            OptimParapm.NN = NN;
            
            OptimParapp.nq = obj.nq/4;
            OptimParamp.nq = obj.nq/4;
            OptimParamm.nq = obj.nq/4;
            OptimParapm.nq = obj.nq/4;
            
        end
        
        function [psi,out] = ReconstructDistribution(obj,u)
            persistent betapp BMpp betamp BMmp betamm BMmm betapm BMpm
            
            if isempty(betapp) || ~isequal(size(betapp),[size(u,1)/4,size(u,2)])
                betapp = repmat([log(1/pi);zeros(size(u,1)/4-1,1)],1,size(u,2));
                betamp = betapp;
                betamm = betapp;
                betapm = betapp;
                BMpp = cell(size(u,2),1);
                [BMpp{:}] = deal(eye(size(u,1)/4));
                BMmp = BMpp;
                BMmm = BMpp;
                BMpm = BMpp;
            end
            
            NQ = obj.nq/4;
            NN = obj.NumberMoments/4;
            OptimParapp = obj.OptimParapp; %Throw out what is not necessary
            OptimParapp.w = OptimParapp.w';
            
            OptimParamp = obj.OptimParamp; %Throw out what is not necessary
            OptimParamp.w = OptimParamp.w';
            
            OptimParamm = obj.OptimParamm; %Throw out what is not necessary
            OptimParamm.w = OptimParamm.w';
            
            OptimParapm = obj.OptimParapm; %Throw out what is not necessary
            OptimParapm.w = OptimParapm.w';
            
            psi = zeros(size(obj.BasisFunctionsAtQuadrature,3),size(u,2));
            for i=1:size(u,2);
                Upp = u(1:NN,i);
                Ump = u(NN+1:2*NN,i);
                Umm = u(2*NN+1:3*NN,i);
                Upm = u(3*NN+1:4*NN,i);
                [alphaF,betapp(:,i),BMpp{i}] = dualAdaptivePoly(Upp, betapp(:,i) ,BMpp{i},OptimParapp);
                psi(1:NQ,i) = Upp(1)*exp(alphaF'*OptimParapp.p);
                [alphaF,betamp(:,i),BMmp{i}] = dualAdaptivePoly(Ump, betamp(:,i) ,BMmp{i},OptimParamp);
                psi(NQ+1:2*NQ,i) = Ump(1)*exp(alphaF'*OptimParamp.p);
                [alphaF,betamm(:,i),BMmm{i}] = dualAdaptivePoly(Umm, betamm(:,i) ,BMmm{i},OptimParamm);
                psi(2*NQ+1:3*NQ,i) = Umm(1)*exp(alphaF'*OptimParamm.p);
                [alphaF,betapm(:,i),BMpm{i}] = dualAdaptivePoly(Upm, betapm(:,i) ,BMpm{i},OptimParapm);
                psi(3*NQ+1:4*NQ,i) = Upm(1)*exp(alphaF'*OptimParapm.p);
            end
            out = [];
            
            
        end
        
        function [UC,J] = NeutronKernel(obj,~,U)
            NN = obj.NumberMoments/4;
            
            UC = bsxfun(@times,obj.phiIso/obj.Density(obj.phiIso),obj.Density(U))-U;
            blk1 = [obj.phiIso/obj.Density(obj.phiIso),zeros(obj.NumberMoments,NN-1)];
            J1 = [blk1,blk1,blk1,blk1];
            J = repmat(J1-diag(ones(obj.NumberMoments,1)),[1 1 size(UC,2)]);
        end
        
        
        function [Fxp,Fxm,Fyp,Fym,Uc,Jc] = Closure(obj,U,R,scatteringtrue)
            psi = obj.ReconstructDistribution(reshape(U,size(U,1),[]));
            if scatteringtrue
                [Uc,Jc] = R.Collision(psi,U);
                Uc = reshape(Uc,size(U));
            else
                Uc = 0*U;
                Jc = zeros(obj.NumberMoments,obj.NumberMoments,size(U,2));
            end
            
            psixp = bsxfun(@times,psi,[obj.Omegaxpp;0*obj.Omegaxmp;0*obj.Omegaxmm;obj.Omegaxpm]);
            psiyp = bsxfun(@times,psi,[obj.Omegaypp;obj.Omegaymp;0*obj.Omegaymm;0*obj.Omegaypm]);
            psixm = bsxfun(@times,psi,[0*obj.Omegaxpp;obj.Omegaxmp;obj.Omegaxmm;0*obj.Omegaxpm]);
            psiym = bsxfun(@times,psi,[0*obj.Omegaypp;0*obj.Omegaymp;obj.Omegaymm;obj.Omegaypm]);
            Fxp = obj.ProjectToBasis(psixp);
            Fxp = reshape(Fxp,size(U));
            Fxm = obj.ProjectToBasis(psixm);
            Fxm = reshape(Fxm,size(U));
            Fyp = obj.ProjectToBasis(psiyp);
            Fyp = reshape(Fyp,size(U));
            Fym = obj.ProjectToBasis(psiym);
            Fym = reshape(Fym,size(U));
        end
        
        function LB = LaplaceBeltrami(obj,psi)
            error('not implemented for quarter moments')
        end
        
        function u = ProjectToBasis(obj,psi)
            I=1;
            NN = obj.NumberMoments/4;
            nq = obj.nq/4;
            if ~(size(psi,1)==obj.nq)
                I = find(size(psi)==obj.nq);
                psi = permute(psi,[I,setdiff(1:ndims(psi),I)]);
            end
            sizenew = size(psi);
            upp = projectBasis(psi(1:nq,:),obj.QuadratureWeights_pp,obj.BasisFunctionsAtQuadratureSqueezed(1:NN,:),3);
            ump = projectBasis(psi(nq+1:2*nq,:),obj.QuadratureWeights_mp,obj.BasisFunctionsAtQuadratureSqueezed(NN+1:2*NN,:),3);
            umm = projectBasis(psi(2*nq+1:3*nq,:),obj.QuadratureWeights_mm,obj.BasisFunctionsAtQuadratureSqueezed(2*NN+1:3*NN,:),3);
            upm = projectBasis(psi(3*nq+1:4*nq,:),obj.QuadratureWeights_pm,obj.BasisFunctionsAtQuadratureSqueezed(3*NN+1:4*NN,:),3);
            u = reshape([upp;ump;umm;upm],[4*NN,sizenew(2:end)]);
            u = permute(u,[2:I,1,I+1:ndims(u)]);
        end
        
        function ModelData = save(obj,varargin)
            warning('off','MATLAB:structOnObject');
            ModelData = struct(obj);
            ModelData.Name = class(obj);
            warning('on','MATLAB:structOnObject');
            if nargin==3
                path = varargin{1};
                filename = varargin{2};
                save([path filename],'ModelData')
            end
        end
        
        function obj = load(obj,varargin)
            if nargin==2
                ModelData = varargin{1};
            else
                path = varargin{1};
                filename = varargin{2};
                load([path filename]); %gives ModelData
            end
            if strcmp(class(obj),ModelData.Name)
                ModelData = rmfield(ModelData,'Name');
                F = fieldnames(ModelData);
                %try
                for i=1:length(F)
                    obj.(F{i}) = ModelData.(F{i});
                end
                %catch
                %     error('It appears that you wanted to load old data')
                % end
            else
                error(sprintf('It appears that you wanted to load the data for class "%s" into class "%s"\n',ModelData.Name,class(obj))); %#ok<SPERR>
            end
            
        end
    end
end

