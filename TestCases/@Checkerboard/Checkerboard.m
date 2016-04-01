classdef Checkerboard < TestCase
    %Checkerboard: Children of the TestCase class.
    %For more information see also TestCase.
    %     
    % (c) Florian Schneider (schneider@mathematik.uni-kl.de) 2016
    % See licence.txt for more information about usage and redistribution
    
    properties
    end
    
    methods
        
        function obj = Checkerboard(varargin)
            
            %Inputparser
            p = inputParser;
            p.KeepUnmatched=true;
            p.parse(varargin{:});
            
            InputStruct = p.Unmatched;
            
            InputStruct.InitialCondition = @(x,y) 1e-10/4/pi+0*x;
            InputStruct.BoundaryCondition_handle = @(x,y,t) 1e-10/4/pi+0*x;
            
            %Absorption
            Spot1 = @(x,y) double(x<=2 & x>= 1 & y<=2 & y>=1);
            Spot2 = @(x,y) double(x<=2 & x>= 1 & y<=4 & y>=3);
            Spot3 = @(x,y) double(x<=2 & x>= 1 & y<=6 & y>=5);
            Spot4 = @(x,y) double(x<=6 & x>= 5 & y<=2 & y>=1);
            Spot5 = @(x,y) double(x<=6 & x>= 5 & y<=4 & y>=3);
            Spot6 = @(x,y) double(x<=6 & x>= 5 & y<=6 & y>=5);
            Spot7 = @(x,y) double(x<=3 & x>= 2 & y<=3 & y>=2);
            Spot8 = @(x,y) double(x<=3 & x>= 2 & y<=5 & y>=4);
            Spot9 = @(x,y) double(x<=5 & x>= 4 & y<=3 & y>=2);
            Spot10 = @(x,y) double(x<=5 & x>= 4 & y<=5 & y>=4);
            Spot11 = @(x,y) double(x<=4 & x>= 3 & y<=2 & y>=1);
            
            Spot0 = @(x,y) double(y<=4 & y>= 3 & x<=4 & x>= 3);
            InputStruct.sigma_a = @(x,y,t) 10*(Spot1(x,y)+Spot2(x,y)+Spot3(x,y)+Spot4(x,y)+Spot5(x,y)+Spot6(x,y)+Spot7(x,y)+Spot8(x,y)+Spot9(x,y)+Spot10(x,y)+Spot11(x,y));
            %Scattering
            InputStruct.sigma_s = @(x,y,t) 1-(Spot1(x,y)+Spot2(x,y)+Spot3(x,y)+Spot4(x,y)+Spot5(x,y)+Spot6(x,y)+Spot7(x,y)+Spot8(x,y)+Spot9(x,y)+Spot10(x,y)+Spot11(x,y));
            InputStruct.Q = @(x,y,t) Spot0(x,y)*1/4/pi;
            
            
            
            
            obj = obj@TestCase(InputStruct);
            
            
        end
        
        function G = generateGrid(obj,nx,ny,order,~,MinModConstant)
            %% Grid Generator
            
            Params.box = [0,7,0,7];
            
            Params.nx = nx;
            Params.ny = ny;
            Params.order=order;
            Params.MinModConstant = MinModConstant;
            %             Params.BoundaryMode = 'Dirichlet';
            G=Grid(Params);
        end
        
    end
    
    
    
    
    
end

