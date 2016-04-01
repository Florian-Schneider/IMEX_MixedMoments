classdef TwoBeamsOpposing < TwoBeams
    %TwoBeamsOpposing: Children of the TestCase class.
    %For more information see also TestCase.
    %     
    % (c) Florian Schneider (schneider@mathematik.uni-kl.de) 2016
    % See licence.txt for more information about usage and redistribution
    
    properties
    end
    
    methods
        
        function obj = TwoBeamsOpposing(varargin)
            
            %Inputparser
            p = inputParser;
            p.KeepUnmatched=true;
            p.parse(varargin{:});
            
            InputStruct = p.Unmatched;

            obj = obj@TwoBeams(InputStruct);
            
            obj.Spot1 = @(x,y) double(y>=0.45 & y<= 0.55 & x<=0);
            obj.Spot2 = @(x,y) double(y>=0.45 & y<= 0.55 & x>=1);
            obj.mu1 = 0;
            obj.mu2 = 0;
            obj.phi1 = 0;
            obj.phi2 = pi;
            
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

