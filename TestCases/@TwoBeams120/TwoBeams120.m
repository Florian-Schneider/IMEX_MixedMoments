classdef TwoBeams120 < TestCase
    %TwoBeams120: Children of the TestCase class.
    %For more information see also TestCase.
    %     
    % (c) Florian Schneider (schneider@mathematik.uni-kl.de) 2016
    % See licence.txt for more information about usage and redistribution
    
    properties
    end
    
    methods
        
        function obj = TwoBeams120(varargin)
            
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
            
            
            
            Psi0 = 100/4/pi;
            sigma = sqrt(0.05); %Width of the beam
            
            InputStruct.BoundaryCondition_handle = @(mu,phi,x,y,t) Psi0*(double(x<=0.05 & y<=0.05).*(-1e-4/4/pi+exp(-(mu.^2/2/sigma^2)-((phi-pi/4).^2/2/sigma^2))+exp(-(mu.^2/2/sigma^2)-((phi-(2.*pi+pi/4)).^2/2/sigma^2)))...
                                                         +double(y>=0.45 & y<= 0.55 & x>=1).*(-1e-4/4/pi+exp(-(mu.^2/2/sigma^2)-((phi-pi).^2/2/sigma^2))+exp(-(mu.^2/2/sigma^2)-((phi-(2.*pi+pi)).^2/2/sigma^2))))...
                                                         +1e-4/4/pi;
            %
            InputStruct.BoundaryCondition_constant = true;
            obj = obj@TestCase(InputStruct);
            
            
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
    end
    
    
    
    
    
end

