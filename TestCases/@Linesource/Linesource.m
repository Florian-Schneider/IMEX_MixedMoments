classdef Linesource < TestCase
    %Linesource: Children of the TestCase class.
    %For more information see also TestCase.
    %     
    % (c) Florian Schneider (schneider@mathematik.uni-kl.de) 2016
    % See licence.txt for more information about usage and redistribution
    
    properties
    end
    
    methods
        
        function obj = Linesource(varargin)
            
            %Inputparser
            p = inputParser;
            p.KeepUnmatched=true;
            p.parse(varargin{:});         
            
            InputStruct = p.Unmatched;
            if ~isfield(p.Unmatched,'InitialCondition')
                sigma = 0.03;
%                 InputStruct.InitialCondition = @(mu,phi,x,y) 1/4/pi*max(exp(-10*((x).^2+y.^2)/sigma^2),1e-4);
                InputStruct.InitialCondition = @(x,y) max(1/(8*pi*sigma^2)*exp(-(x.^2+y.^2)/2/sigma^2),1e-4); %Code identifies the arguments (x,y) to calculate isotropic conditions
            end
            if ~isfield(p.Unmatched,'sigma_a')
                InputStruct.sigma_a = @(x,y,t) 0*x;
            end
            if ~isfield(p.Unmatched,'sigma_s')
                InputStruct.sigma_s = @(x,y,t) 1+0*x;
            end
            if ~isfield(p.Unmatched,'Q')
                InputStruct.Q = @(x,y,t) 0*x;
            end
            
            %InputStruct.BoundaryCondition_handle = @(mu,phi,x,y,t) 1/4/pi*max(exp(-10*((x).^2+y.^2)/sigma^2),1e-4);
            InputStruct.BoundaryCondition_handle = @(x,y,t) 1/4/pi*max(exp(-10*((x).^2+y.^2)/sigma^2),1e-4);
            
            obj = obj@TestCase(InputStruct);
            
            
        end
        
        function G = generateGrid(obj,nx,ny,order,~,MinModConstant)
            %% Grid Generator
            
            Params.box = [-0.5,0.5,-0.5,0.5];
            
            Params.nx = nx;
            Params.ny = ny;
            Params.order=order;
            Params.MinModConstant = MinModConstant;
%             Params.BoundaryMode = 'Dirichlet';
            G=Grid(Params);    
        end
        
    end
    
    
    
        
    
end

