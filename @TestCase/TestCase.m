classdef TestCase < handle
    %TestCase: TestCase class for the IMEX scheme handling the test case.
    %Consists of the following information:
    % t_final: end time of the simulation
    % sigma_a: handle containing the spatial dependence on the absorption
    %    coefficient
    % sigma_s: handle containing the spatial dependence on the scattering
    %    coefficient
    % Q: handle containing the spatial and angular dependence on the source
    % InitialCondition: handle containing the spatial and angular dependence on the
    % initial condition for the underlying kinetic distribution psi
    %
    % For examples see also the Tests folder
    % 
    % (c) Florian Schneider (schneider@mathematik.uni-kl.de) 2016
    % See licence.txt for more information about usage and redistribution
    
    properties
        
        sigma_a_handle;
        sigma_s_handle;
        Q_handle;
        sigma_a;
        sigma_s;
        Q;
        sigma_a_constant;
        sigma_s_constant;
        Q_constant;
        t_final;
        
        InitialCondition;
        BoundaryCondition_handle;
        BoundaryCondition;
        BoundaryCondition_constant; %true if the b.c. are constant over time
        
    end
    
    methods
        
        function obj = TestCase(varargin)
            
            %Inputparser
            p = inputParser;
            p.addParameter('t_final',1);
            p.addParameter('sigma_a',@(x,y) 0*x);
            p.addParameter('sigma_s',@(x,y) 0*x);
            p.addParameter('Q',@(x,y,mu,phi) 0*x);
            p.addParameter('InitialCondition',@(mu,phi,x,y) 1+0*x);
            p.addParameter('BoundaryCondition_handle',@(mu,phi,x,y,t) 0*x);
            p.addParameter('BoundaryCondition_constant',true);
            p.addParameter('sigma_a_constant',true);
            p.addParameter('sigma_s_constant',true);
            p.addParameter('Q_constant',true);
            p.addParameter('Path',[]);
            p.addParameter('Filename',[]);
            p.addParameter('TestCaseData',[]);
            
            p.parse(varargin{:});
            
            if (isempty(p.Results.Filename) || not(exist([p.Results.Path p.Results.Filename],'file')==2)) && isempty(p.Results.TestCaseData)
                
                obj.t_final=p.Results.t_final;
                obj.sigma_a_handle=p.Results.sigma_a;
                obj.Q_handle=p.Results.Q;
                obj.sigma_s_handle=p.Results.sigma_s;
                obj.InitialCondition=p.Results.InitialCondition;
                obj.BoundaryCondition_handle=p.Results.BoundaryCondition_handle;
                obj.BoundaryCondition_constant=p.Results.BoundaryCondition_constant;
                obj.BoundaryCondition = [];
                obj.sigma_a_constant=p.Results.sigma_a_constant;
                obj.sigma_s_constant=p.Results.sigma_s_constant;
                obj.Q_constant=p.Results.Q_constant;
            else
                if isempty(p.Results.TestCaseData)
                    obj.load(p.Results.Path,p.Results.Filename);
                else
                    obj.load(p.Results.TestCaseData);
                end
            end
        end
        
        function G = generateGrid(obj,dx,dy,order,MomentModel,MinModConstant)
            G = dx*dy; %#ok<NASGU>
            error('TestCase:generateGrid_Dummy','generateGrid has to be created in every single TestCase Subclass!');
            
            
        end
        
        function [BoundaryCondition,BoundaryCondition_psi] = getBC2(obj,G,M,t) %Grid and MomentModel Input
            if (~obj.BoundaryCondition_constant || (t==0 && isempty(obj.BoundaryCondition))) && ~isempty(obj.BoundaryCondition_handle)
                BoundaryCondition_psi = zeros(M.nq,G.nx(3));
                mu = M.QuadraturePoints(1,:);
                phi = M.QuadraturePoints(2,:);
                for i=1:M.nq
                    handle = @(x,y) obj.BoundaryCondition_handle(mu(i),phi(i),x,y,t);
                    BoundaryCondition_psi(i,:) = G.compute_cell_averages(handle,1:G.nx(3));
                end
                BoundaryCondition = M.ProjectToBasis(BoundaryCondition_psi);
            else
                BoundaryCondition = obj.BoundaryCondition;
            end
            
            obj.BoundaryCondition = BoundaryCondition;
            
        end
        
        function [BoundaryCondition] = getBC(obj,G,M,t) %Grid and MomentModel Input
            if (~obj.BoundaryCondition_constant || t==0 ) && ~isempty(obj.BoundaryCondition_handle)
                
                if strcmp(G.BoundaryMode,'periodic')
                    psi = zeros(length(G.BoundaryCells),size(G.Basis,3),M.nq);
                    psival = sum(bsxfun(@times,permute(psi,[1 3 4 5 2]),permute(G.Basis,[4 5 1 2 3])),5); 
                    BoundaryConditionVal = M.ProjectToBasis(psival);
                else
                
                if isIsotropic(obj.BoundaryCondition_handle)
                    psi = repmat(G.projectToBasis(@(x,y) obj.BoundaryCondition_handle(x,y,t),G.BoundaryCells),1,1,M.nq);
                else
                    psi = zeros(length(G.BoundaryCells),size(G.Basis,3),M.nq);
                    for i=1:M.nq
                        psi(:,:,i) = G.projectToBasis(@(x,y) obj.BoundaryCondition_handle(M.QuadraturePoints(1,i),M.QuadraturePoints(2,i),x,y,t),G.BoundaryCells); %Very slow...
                    end
                end
                
                psival = sum(bsxfun(@times,permute(psi,[1 3 4 5 2]),permute(G.Basis,[4 5 1 2 3])),5); %TODO warum ist psival nicht konstant entlang der QP?
                psibar = permute(psi(:,1,:),[1 3 2]);
                theta = (bsxfun(@times,psibar,1./bsxfun(@plus,psibar,-psival)));
                theta(theta>1 | theta < 0) = 1;
                theta(theta<1) = max(theta(theta<1),0);
                theta2 = min(min(theta,[],4),[],3);
                psival2 = bsxfun(@plus,(1-theta2).*psibar,bsxfun(@times,theta2,psival));
                BoundaryConditionVal = M.ProjectToBasis(psival2);
                end
                
                nBasis = size(G.Basis,3);
                s = size(psi);
                s(3) = M.NumberMoments;
                BoundaryCondition = zeros(s);
                for i=1:nBasis
                    BoundaryCondition(:,i,:) = (G.massmatrix(i)).^(-1)*sum(sum(bsxfun(@times,BoundaryConditionVal,permute(G.Quadrature.W.*G.Basis(:,:,i),[3 4 1 2])),4),3);
                end
            else
                BoundaryCondition = obj.BoundaryCondition;
            end
            
            obj.BoundaryCondition = BoundaryCondition;
            
        end
        
        function [sigma_s] = getSigma_s(obj,G,t) %Grid and MomentModel Input
            if (~obj.sigma_s_constant || (t==0 || isempty(obj.sigma_s))) && ~isempty(obj.sigma_s_handle)
                sigma_s=G.projectToBasis(@(x,y) obj.sigma_s_handle(x,y,t),G.InteriorCells);
                sigma_s(:,2:end) = 0; %only cell-means
            else
                sigma_s = obj.sigma_s;
            end
            obj.sigma_s = sigma_s;
        end
        
        function [sigma_a] = getSigma_a(obj,G,t) %Grid and MomentModel Input
            if (~obj.sigma_a_constant || (t==0 || isempty(obj.sigma_a))) && ~isempty(obj.sigma_a_handle)
                sigma_a=G.projectToBasis(@(x,y) obj.sigma_a_handle(x,y,t),G.InteriorCells);
                sigma_a(:,2:end) = 0; %only cell-means
            else
                sigma_a = obj.sigma_a;
            end
            obj.sigma_a = sigma_a;
        end
        
        function [Q] = getQ(obj,G,M,t) %Grid and MomentModel Input
            if (~obj.Q_constant || (t==0 || isempty(obj.Q))) && ~isempty(obj.Q_handle)
                if isIsotropic(obj.Q_handle)
                    psi = repmat(G.projectToBasis(@(x,y) obj.Q_handle(x,y,t),G.InteriorCells),1,1,M.nq);
                else
                    psi = zeros(length(G.InteriorCells),size(G.Basis,3),M.nq);
                    for i=1:M.nq
                        psi(:,:,i) = G.projectToBasis(@(x,y) obj.Q_handle(M.QuadraturePoints(1,i),M.QuadraturePoints(2,i),x,y,t),G.InteriorCells); %Very slow...
                    end
                end
                Q = M.ProjectToBasis(psi);
            else
                Q = obj.Q;
            end
            
            obj.Q = Q;
        end
        
        
        function TestCaseData = save(obj,varargin)
            warning('off','MATLAB:structOnObject');
            TestCaseData = struct(obj);
            TestCaseData.Name = class(obj);
            warning('on','MATLAB:structOnObject');
            if nargin==3
                path = varargin{1}; 
                filename = varargin{2}; 
                save([path filename],'TestCaseData')
            end
        end
        
        function load(obj,varargin)
            if nargin==2
                TestCaseData = varargin{1}; 
            else
                path = varargin{1}; 
                filename = varargin{2}; 
                load([path filename]); %gives TestCaseData
            end
            if strcmp(class(obj),TestCaseData.Name)
                TestCaseData = rmfield(TestCaseData,'Name');
                F = fieldnames(TestCaseData);
                %try
                for i=1:length(F)
                    obj.(F{i}) = TestCaseData.(F{i});
                end
                %catch 
               %     error('It appears that you wanted to load old data')
               % end
            else
                error(sprintf('It appears that you wanted to load the data for class "%s" into class "%s"\n',TestCaseData.Name,class(obj))); %#ok<SPERR>
            end
            
        end
    end
    
    
    
end
