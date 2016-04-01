classdef MomentModel < handle
    %MomentModel: Super class for all moment models. Consists of 
    %  the MomentOrder N and (if necessary) the quadrature order to
    %  approximate spherical integrals. Only invoked by the children
    %  classes.
    %
    % (c) Florian Schneider (schneider@mathematik.uni-kl.de) 2016
    % See licence.txt for more information about usage and redistribution
    %
    
    
    properties
        BasisFunctions
        BasisFunctionsAtQuadrature
        BasisFunctionsAtQuadratureSqueezed
        MomentOrder
        NumberMoments
        nq
        OptimPara
        Omegax
        Omegay
        phiIso
        positivityFlag
        QuadraturePoints
        QuadratureWeights
        QuadratureWeightsSqueezed
        Density
        
        TimeStepping
        
        ObviouslyRealizableCounter
        
        QuadratureOrder
    end
    
    methods
        
        
        function obj = MomentModel(varargin)
            %Inputparser
            p = inputParser;
            p.addParameter('MomentOrder',[]);
             p.addParameter('QuadratureOrder',20);
            
            p.addParameter('QuadraturePoints',[]);
            p.addParameter('QuadratureWeights',[]);
            
            p.addParameter('Path',[]);
            p.addParameter('Filename',[]);
            p.addParameter('ModelData',[]);
            
            p.parse(varargin{:});
            
            if (isempty(p.Results.Filename) || not(exist([p.Results.Path p.Results.Filename],'file')==2)) && isempty(p.Results.ModelData)
            
            
            obj.QuadratureOrder = p.Results.QuadratureOrder;
            
            if ~isempty(p.Results.QuadraturePoints)
                obj.QuadraturePoints=p.Results.QuadraturePoints;
                obj.QuadratureWeights=reshape(p.Results.QuadratureWeights,1,1,[]); %1x1xnq
                obj.QuadratureWeightsSqueezed=p.Results.QuadratureWeights;
            else
                Q = LebedevQuadrature(obj.QuadratureOrder);
                obj.QuadraturePoints = Q(1:2,:);
                w = Q(3,:);
                obj.QuadratureWeights=reshape(w,1,1,[]); %1x1xnq
                obj.QuadratureWeightsSqueezed=w;
            end
            obj.nq = size(obj.QuadraturePoints,2); %QuadraturePoints = 2xnq
                
            obj.MomentOrder=p.Results.MomentOrder;
            [ obj.BasisFunctions,obj.BasisFunctionsAtQuadrature,obj.BasisFunctionsAtQuadratureSqueezed] = obj.generateBasisFunctions; %@(mu,phi) [BF1(mu,phi);BF2(mu,phi);...]
            
            obj.NumberMoments = size(obj.BasisFunctionsAtQuadrature,1);
            
            obj.Omegax = transpose(sqrt(1-obj.QuadraturePoints(1,:).^2).*cos(obj.QuadraturePoints(2,:)));
            obj.Omegay = transpose(sqrt(1-obj.QuadraturePoints(1,:).^2).*sin(obj.QuadraturePoints(2,:)));
            
            
            psi = ones(obj.nq,1);
            obj.phiIso = obj.ProjectToBasis(psi);
            obj.OptimPara = generateOptimPara(obj);
            obj.positivityFlag = 1;
            sizecap = @(s) s(2:end);
            obj.Density = @(U) reshape(U(1,:),[1,sizecap(size(U))]);
            else
                if isempty(p.Results.ModelData)
                    obj = obj.load(p.Results.Path,p.Results.Filename);
                else
                    obj = obj.load(p.Results.ModelData);
                end
            end
            obj.TimeStepping.A = 1;
            obj.TimeStepping.b = 1;
            obj.TimeStepping.Stage = 1;
            obj.TimeStepping.dt = 0;
        end
        
                
        function a = ReconstructDistribution(obj,u)%#ok<INUSD,STOUT>
            error('Abstract class, must be implemented in subclass!');
            a = [];
            
        end
        
        function BF = generateBasisFunctions(obj)%#ok<MANU,STOUT>
            error('Abstract class, must be implemented in subclass!');
            BF = [];
        end
        
        function [LB,J] = LaplaceBeltrami(obj,psi) %#ok<INUSD,STOUT>
            error('Abstract class, must be implemented in subclass!');
        end
        
        function OptimPara = generateOptimPara(obj)
            OptimPara.itermax = 0;
            OptimPara.rList = 0;%[0,1e-4,1e-2,1e-1];
            OptimPara.tol = 1e-6; % tolerance on the norm of the gradient in the
            %                       stopping criterion
            OptimPara.rList = [0, 1e-8, 1e-6, 1e-4,1e-2,1e-1];% series of regularization values to try
            OptimPara.k0 = 50; % number of Newton iterations to try before increasing r
            OptimPara.itermax = max(OptimPara.itermax, OptimPara.k0 * length(OptimPara.rList)); % maximum number of iterations
            OptimPara.rMax = max(OptimPara.rList);
%             r0 = 1e-4; % should be positive and fairly small
%             kappa = 0.5; % should be in ]0, 1[
            r0 = 1e-6; % should be positive and fairly small
            kappa = 0.1; % should be in ]0, 1[

            ellMax = ceil(log(eps / r0) / log(kappa));
            OptimPara.rList2 = [r0 * kappa.^(0:ellMax), 0];
            OptimPara.p = obj.BasisFunctionsAtQuadratureSqueezed;
            OptimPara.phiIso = obj.phiIso;
            OptimPara.w = obj.QuadratureWeightsSqueezed';
            OptimPara.NN = obj.NumberMoments;
            OptimPara.nq = size(obj.BasisFunctionsAtQuadratureSqueezed,2);
            
            
        end

        function [Fxp,Fxm,Fyp,Fym,Uc,Jc] = Closure(obj,U,R,scatteringtrue)
            [psi,out] = obj.ReconstructDistribution(reshape(U,size(U,1),[]));
            if ~isempty(out)
                if ~isfield(R.Statistics,'Optimization')
                    R.Statistics.Optimization = repmat(out,[1,length(R.t_Frame),length(obj.TimeStepping.A)]);
                end
                R.Statistics.Optimization(:,R.Frame_cnt,obj.TimeStepping.Stage) = out;
            else
                R.Statistics.Optimization = [];
            end
            if scatteringtrue
                [Uc,Jc] = R.Collision(psi,U);
                Uc = reshape(Uc,size(U));
            else
                Uc = 0*U;
                Jc = zeros(obj.NumberMoments,obj.NumberMoments,size(U,2));
            end
            
            psixp = bsxfun(@times,psi,obj.Omegax.*(obj.Omegax>=0));
            psixm = bsxfun(@times,psi,obj.Omegax.*(obj.Omegax<=0));
            psiyp = bsxfun(@times,psi,obj.Omegay.*(obj.Omegay>=0));
            psiym = bsxfun(@times,psi,obj.Omegay.*(obj.Omegay<=0));
            Fxp = obj.ProjectToBasis(psixp);
            Fxp = reshape(Fxp,size(U));
            Fxm = obj.ProjectToBasis(psixm);
            Fxm = reshape(Fxm,size(U));
            Fyp = obj.ProjectToBasis(psiyp);
            Fyp = reshape(Fyp,size(U));
            Fym = obj.ProjectToBasis(psiym);
            Fym = reshape(Fym,size(U));                        
        end
        
        function u = ProjectToBasis(obj,psi)
            I=1;
            if ~(size(psi,1)==obj.nq)
                I = find(size(psi)==obj.nq);
                psi = permute(psi,[I,setdiff(1:ndims(psi),I)]);
            end
            sizenew = size(psi);
            u = projectBasis(psi,obj.QuadratureWeightsSqueezed,obj.BasisFunctionsAtQuadratureSqueezed,3);
            u = reshape(u,[size(u,1),sizenew(2:end)]);
            u = permute(u,[2:I,1,I+1:ndims(u)]);
%             if ~isequal(size(psi,2),obj.nq)
%                 if isequal(size(psi,1),obj.nq)
%                     psi = permute(psi,[3 2 1]);
%                 else
%                     error('Argument 2 (psi) does not have the right dimension')
%                 end
%             else
%                 psi = permute(psi,[3 1 2]);
%             end
%             
%             u = sum(bsxfun(@times,permute(obj.QuadratureWeights,[1 3 2]),bsxfun(@times,obj.BasisFunctionsAtQuadrature,psi)),3);

%             u = sum(bsxfun(@times,bsxfun(@times,obj.QuadratureWeights,obj.BasisFunctionsAtQuadrature),permute(psi,[3 2 1])),3);
        end
        
        function Uout = SlopeLimiter(~,U,G)
           tol = 1e-5;
           Uout = U(G.InteriorCells,:,:);
           for i=1:size(U,3)
               Ubar = U(G.InteriorCells,1,i);
               Ux = U(G.InteriorCells,G.BasisXdiff,i);
               Uy = U(G.InteriorCells,G.BasisYdiff,i);
               Uleft = U(G.LeftNeighbor,1,i);
               Uright = U(G.RightNeighbor,1,i);
               Ulower = U(G.LowerNeighbor,1,i);
               Uupper = U(G.UpperNeighbor,1,i);
               Uxmod = TVBminmod(Ux,Uright-Ubar,Ubar-Uleft,G.dx,G.MinModConstant);
               Uymod = TVBminmod(Uy,Uupper-Ubar,Ubar-Ulower,G.dy,G.MinModConstant);
               I = abs(Ux-Uxmod)>tol | abs(Uy-Uymod)>tol;
               Uout(:,G.BasisXdiff,i) = Uxmod;
               Uout(:,G.BasisYdiff,i) = Uymod;
               Uout(I,setdiff(1:size(G.Basis,3),[1,G.BasisXdiff,G.BasisYdiff]),i) = 0;
           end
        end
        
        function [U,RealizabilityCounter] = RealizabilityLimiter(obj,U,Basis)
            RealizabilityCounter = [0 0];
        end
        
        function u0 = ProjectToDensity(obj,U,c)
            %Moments are stored in the c-th component of U. Function
            %projects to u0
            Up = permute(U,[c,1:c-1,c+1:length(size(U))]);
            u0 = Up(1,:)'/obj.BasisFunctionsAtQuadratureSqueezed(1);
            s = size(Up);
            u0 = reshape(u0,s(2:end));
            
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
