classdef TwoBeams < TestCase
    %TwoBeams: Children of the TestCase class.
    %For more information see also TestCase.
    %
    % (c) Florian Schneider (schneider@mathematik.uni-kl.de) 2016
    % See licence.txt for more information about usage and redistribution
    
    properties
        Spot1
        Spot2
        mu1
        mu2
        phi1
        phi2
    end
    
    methods
        
        function obj = TwoBeams(varargin)
            
            %Inputparser
            p = inputParser;
            p.KeepUnmatched=true;
            p.parse(varargin{:});
            
            InputStruct = p.Unmatched;
            if ~isfield(p.Unmatched,'InitialCondition')
                InputStruct.InitialCondition = @(x,y) 1e-4/4/pi;
            end
            if ~isfield(p.Unmatched,'sigma_a')
                InputStruct.sigma_a = @(x,y,t) 0*x;
            end
            if ~isfield(p.Unmatched,'sigma_s')
                InputStruct.sigma_s = @(x,y,t) 0*x;
            end
            if ~isfield(p.Unmatched,'Q')
                InputStruct.Q = @(x,y,t) 0*x;
            end
            InputStruct.BoundaryCondition_handle = @(mu,phi,x,y,t) NaN;
            %
            InputStruct.BoundaryCondition_constant = true;
            obj = obj@TestCase(InputStruct);
            
            obj.Spot1 = @(x,y) double(y>=0.45 & y<= 0.55 & x<=0);
            obj.Spot2 = @(x,y) double(x>=0.45 & x<= 0.55 & y<=0);
            obj.mu1 = 0;
            obj.mu2 = 0;
            obj.phi1 = 0;
            obj.phi2 = pi/2;
            
            
        end
        
        function G =  generateGrid(obj,nx,ny,order,~,MinModConstant)
            %% Grid Generator
            
            Params.box = [0,1,0,1];
            
            Params.nx = nx;
            Params.ny = ny;
            Params.order=order;
            Params.MinModConstant = MinModConstant;
            %             Params.BoundaryMode = 'Dirichlet';
            G=Grid(Params);
        end
        
        
        function [BoundaryCondition] = getBC(obj,G,M,t) %Grid and MomentModel Input
            if (~obj.BoundaryCondition_constant || t==0 ) && ~isempty(obj.BoundaryCondition_handle)
                
                S1 = G.projectToBasis(obj.Spot1,G.BoundaryCells);
                S1 = S1>0;
                S2 = G.projectToBasis(obj.Spot2,G.BoundaryCells);
                S2 = S2>0;
                
                if isa(M,'QK1')
                    U1 = [1;1;0;0;0;0;0;0;0;1;1;0];
                    U2 = [1;0;1;1;0;1;0;0;0;0;0;0];
                    BoundaryCondition = bsxfun(@times,S1,permute(U1,[2 3 1]))+bsxfun(@times,S2,permute(U2,[2 3 1]));  
                else
                    eval(sprintf('Mtmp = %s(''MomentOrder'',%d);',class(M),M.MomentOrder));
                    eval(sprintf('Mtmp2 = %s(''MomentOrder'',%d);',class(M),M.MomentOrder));
                    Mtmp.QuadraturePoints = [obj.mu1;obj.phi1];
                    [ Mtmp.BasisFunctions,Mtmp.BasisFunctionsAtQuadrature,Mtmp.BasisFunctionsAtQuadratureSqueezed] = Mtmp.generateBasisFunctions;
                    Mtmp2.QuadraturePoints = [obj.mu2;obj.phi2];
                    [ Mtmp2.BasisFunctions,Mtmp2.BasisFunctionsAtQuadrature,Mtmp2.BasisFunctionsAtQuadratureSqueezed] = Mtmp2.generateBasisFunctions;
                    
                    BoundaryCondition = bsxfun(@times,S1,permute(Mtmp.BasisFunctionsAtQuadrature,[2 3 1]))+bsxfun(@times,S2,permute(Mtmp2.BasisFunctionsAtQuadrature,[2 3 1]));                    
                end
                BoundaryCondition = bsxfun(@plus,BoundaryCondition,1e-4*permute(M.phiIso,[2 3 1])/M.phiIso(1)); 
            else
                BoundaryCondition = obj.BoundaryCondition;
            end
            
            obj.BoundaryCondition = BoundaryCondition;
            
        end
        
    end
    
    
    
    
    
end

