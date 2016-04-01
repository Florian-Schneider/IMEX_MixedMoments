classdef Grid < handle
    %Grid: Grid class for the IMEX scheme handling the underlying mesh and
    %bases. Naturally invoked by the Testcase class.
    % For examples see also the Tests folder, Testcase
    % 
    % (c) Florian Schneider (schneider@mathematik.uni-kl.de) 2016
    % See licence.txt for more information about usage and redistribution
    
    properties
        box                         %[xmin,xmax,ymin,ymax]
        xmin
        xmax
        ymin
        ymax
        dx
        dy
        nx                          %number cells in x direction
        ny                          %number cells in y direction
        nint
        InteriorCells
        BoundaryCells
        LeftNeighbor
        RightNeighbor
        UpperNeighbor
        LowerNeighbor
        
        BasisOrder
        BasisXdiff
        BasisYdiff
        
        w1
        Quadrature                  %Quadrature in reference element
        QuadratureHR                %High order Quadrature in reference element
        Basis                       %Basis polynomials evaluated at the quadrature points
        BasisHR                     %Basis polynomials evaluated at the higher order quadrature points
        BasisGradx                  %Gradient in x direction evaluated at the quadrature points
        BasisGrady                  %Gradient in y direction evaluated at the quadrature points
        BasisLeft
        BasisRight
        BasisLower
        BasisUpper
        order                       %Order of Reconstruction, degree of polynomials +1
        cell_volumes                %volumes
        midpoints
        
        MinModConstant
        
        massmatrix
        
        LeftNormal
        RightNormal
        LowerNormal
        UpperNormal
        
        LeftEdge_Quadrature_Points% Index for the QP on left edge
        RightEdge_Quadrature_Points
        LowerEdge_Quadrature_Points
        UpperEdge_Quadrature_Points
        LeftEdge_Quadrature_Weights
        RightEdge_Quadrature_Weights
        LowerEdge_Quadrature_Weights
        UpperEdge_Quadrature_Weights
        LeftNEdge_Quadrature_Points% Index for the QP on left edges neighbor
        RightNEdge_Quadrature_Points
        LowerNEdge_Quadrature_Points
        UpperNEdge_Quadrature_Points
        LeftNEdge_Quadrature_Weights
        RightNEdge_Quadrature_Weights
        LowerNEdge_Quadrature_Weights
        UpperNEdge_Quadrature_Weights
        
        BoundaryMode                %{'dirichlet','outflow','periodic'}
        PeriodicMapping             %Mapping for periodic data
        
        
        ObviouslyRealizableCounter
    end
    
    methods
        function obj=Grid(varargin)
            p=inputParser;
            p.addParameter('order',1,@(s) s==1);  %Order of reconstruction
            p.addParameter('box',[]);
            p.addParameter('nx',20);
            p.addParameter('ny',20);
            p.addParameter('MinModConstant',30);
            p.addParameter('BoundaryMode','dirichlet'); %WENO
            
            p.addParameter('Path',[]);
            p.addParameter('Filename',[]);
            p.addParameter('GridData',[]);
            
            p.parse(varargin{:});
            
            if (isempty(p.Results.Filename) || not(exist([p.Results.Path p.Results.Filename],'file')==2)) && isempty(p.Results.GridData)
                %constructor
                obj.order=p.Results.order;
                
                obj.MinModConstant = p.Results.MinModConstant;
                
                obj.BoundaryMode = lower(p.Results.BoundaryMode);
                if ~isempty(p.Results.box)
                    obj.box=p.Results.box;
                    obj.xmin = obj.box(1);
                    obj.xmax = obj.box(2);
                    obj.ymin = obj.box(3);
                    obj.ymax = obj.box(4);
                else
                    error('no bounding box specified');
                end
                obj.nx = p.Results.nx;
                obj.ny = p.Results.ny;
                obj.nint = obj.nx*obj.ny;
                obj.InteriorCells = 1:obj.nint;
                obj.dx = (obj.xmax-obj.xmin)/obj.nx;
                obj.dy = (obj.ymax-obj.ymin)/obj.ny;
                [x,y] = meshgrid(obj.xmin + obj.dx*(1/2+(0:(obj.nx-1))),obj.ymin + obj.dy*(1/2+(0:(obj.ny-1))));
                
                obj.midpoints = [x(:),y(:)];
                cell_neighbor_index = zeros(length(obj.midpoints),4); %[x_l,x_r,y_l,y_r]
                
                % Inner points
                cell_neighbor_index(:,1) = (1:obj.nint)' - obj.ny; %left neighbor in x direction
                cell_neighbor_index(:,2) = (1:obj.nint)' + obj.ny; %right neighbor in x direction
                cell_neighbor_index(:,3) = (1:obj.nint)' - 1; %left neighbor in y direction = downwards
                cell_neighbor_index(:,4) = (1:obj.nint)' + 1; %right neighbor in y direction = upwards
                
                % Create boundary:
                facet1 = [x(abs(x-obj.xmin-obj.dx/2)<obj.dx/100)-obj.dx,y(abs(x-obj.xmin-obj.dx/2)<obj.dx/100)];
                facet2 = [x(abs(x-obj.xmax+obj.dx/2)<obj.dx/100)+obj.dx,y(abs(x-obj.xmax+obj.dx/2)<obj.dx/100)];
                facet3 = [x(abs(y-obj.ymin-obj.dy/2)<obj.dy/100),y(abs(y-obj.ymin-obj.dy/2)<obj.dy/100)-obj.dy];
                facet4 = [x(abs(y-obj.ymax+obj.dy/2)<obj.dy/100),y(abs(y-obj.ymax+obj.dy/2)<obj.dy/100)+obj.dy];
                
                % Link facets
                switch obj.BoundaryMode
                    case 'dirichlet'
                        ntmp = obj.nint;
                        cell_neighbor_index((abs(x-obj.xmin-obj.dx/2)<obj.dx/100),1) = ntmp+(1:length(facet1))';
                        ntmp = ntmp + length(facet1);
                        cell_neighbor_index((abs(x-obj.xmax+obj.dx/2)<obj.dx/100),2) = ntmp+(1:length(facet2))';
                        ntmp = ntmp + length(facet2);
                        cell_neighbor_index((abs(y-obj.ymin-obj.dy/2)<obj.dy/100),3) = ntmp+(1:length(facet3))';
                        ntmp = ntmp + length(facet3);
                        cell_neighbor_index((abs(y-obj.ymax+obj.dy/2)<obj.dy/100),4) = ntmp+(1:length(facet4))';
                        ntmp = ntmp + length(facet4);
                        obj.midpoints = [obj.midpoints;facet1;facet2;facet3;facet4];
                        obj.BoundaryCells = obj.nint+1:ntmp;
                    case 'periodic'
                        cell_neighbor_index((abs(x-obj.xmin-obj.dx/2)<obj.dx/100),1) = find((abs(x-obj.xmax+obj.dx/2)<obj.dx/100)); %Left neighbor of left boundary is right boundary
                        cell_neighbor_index((abs(x-obj.xmax+obj.dx/2)<obj.dx/100),2) = find((abs(x-obj.xmin-obj.dx/2)<obj.dx/100)); %right neighbor of right boundary is left boundary
                        cell_neighbor_index((abs(y-obj.ymin-obj.dy/2)<obj.dy/100),3) = find((abs(y-obj.ymax+obj.dy/2)<obj.dy/100)); %Lower...
                        cell_neighbor_index((abs(y-obj.ymax+obj.dy/2)<obj.dy/100),4) = find((abs(y-obj.ymin-obj.dy/2)<obj.dy/100)); %Upper...
                        obj.BoundaryCells = [];
                    otherwise
                        error('not implemented');
                end
                
                obj.LeftNeighbor = cell_neighbor_index(:,1);
                obj.RightNeighbor = cell_neighbor_index(:,2);
                obj.LowerNeighbor = cell_neighbor_index(:,3);
                obj.UpperNeighbor = cell_neighbor_index(:,4);
                
                %Normals
                obj.LeftNormal = [-1,0];
                obj.RightNormal = [1,0];
                obj.LowerNormal = [0,-1];
                obj.UpperNormal = [0,1];
                
                
                %Quadrature
                x = 0;
                w = 1;
                [X,Y] = meshgrid(x);
                [Wx,Wy] = meshgrid(w);
                W = Wx.*Wy;
                obj.Quadrature.X = X;
                obj.Quadrature.Y = Y;
                obj.Quadrature.W = W;
                
                %Basisgeneration
                [I,J] = meshgrid(0:obj.order-1);
                ind = I+J<obj.order;
                I = I(ind); J = J(ind);
                obj.BasisOrder = [I,J];
                p = zeros([size(X),length(I)]);
                for i=1:length(I)
                    [px] = legPoly(2*X,I(i)); %legendre polynomials
                    [py] = legPoly(2*Y,J(i)); %legendre polynomials
                    p(:,:,i) = px.*py;                   
                end
                obj.Basis  = p;
                obj.massmatrix = squeeze(sum(sum(bsxfun(@times,p.^2,W),2),1));
                obj.cell_volumes = obj.dx*obj.dy;
                
                obj.LeftEdge_Quadrature_Points = 1;
                obj.LeftEdge_Quadrature_Weights = obj.dy;
                obj.RightEdge_Quadrature_Points = 1;
                obj.RightEdge_Quadrature_Weights = obj.dy;
                
                obj.LeftNEdge_Quadrature_Points =obj.RightEdge_Quadrature_Points;
                obj.LeftNEdge_Quadrature_Weights = obj.RightEdge_Quadrature_Weights;
                obj.RightNEdge_Quadrature_Points =obj.LeftEdge_Quadrature_Points;
                obj.RightNEdge_Quadrature_Weights = obj.LeftEdge_Quadrature_Weights;
                
                obj.LowerEdge_Quadrature_Points = 1;
                obj.LowerEdge_Quadrature_Weights = obj.dx;
                obj.UpperEdge_Quadrature_Points = 1;
                obj.UpperEdge_Quadrature_Weights = obj.dx;
                
                obj.LowerNEdge_Quadrature_Points =obj.UpperEdge_Quadrature_Points;
                obj.LowerNEdge_Quadrature_Weights = obj.UpperEdge_Quadrature_Weights;
                obj.UpperNEdge_Quadrature_Points =obj.LowerEdge_Quadrature_Points;
                obj.UpperNEdge_Quadrature_Weights = obj.LowerEdge_Quadrature_Weights;
                
                
                p = permute(p,[3,1,2]);
                obj.BasisLeft = p(:,obj.LeftEdge_Quadrature_Points);
                obj.BasisRight = p(:,obj.RightEdge_Quadrature_Points);
                obj.BasisLower = p(:,obj.LowerEdge_Quadrature_Points);
                obj.BasisUpper = p(:,obj.UpperEdge_Quadrature_Points);
                
                
                
                
                %HR quadrature
                
                [x,w] = lgwt(100,-1,1);
                x = x/2; w = w/2;
                [X,Y] = meshgrid(x);
                [Wx,Wy] = meshgrid(w);
                W = Wx.*Wy;
                obj.QuadratureHR.X = X;
                obj.QuadratureHR.Y = Y;
                obj.QuadratureHR.W = W;
                p = zeros([size(X),length(I)]);
                for i=1:length(I)
                    px = legPoly(2*X,I(i)); %legendre polynomials
                    py = legPoly(2*Y,J(i)); %legendre polynomials
                    p(:,:,i) = px.*py;
                end
                obj.BasisHR  = p;                
            else
                if isempty(p.Results.GridData)
                    obj.load(p.Results.Path,p.Results.Filename);
                else
                    obj.load(p.Results.GridData);
                end
            end
        end
        
        function U=projectToBasis(obj,f,cells)
            
            if nargin<3
                cells = 1:size(obj.midpoints,1);
            end
            if isempty(cells)
                U = zeros(1,0);
                return;
            end
            n_cells = length(cells);
            
            Xphys = bsxfun(@plus,obj.midpoints(cells,1),obj.dx*permute(obj.QuadratureHR.X,[3,1,2]));
            Yphys = bsxfun(@plus,obj.midpoints(cells,2),obj.dy*permute(obj.QuadratureHR.Y,[3,1,2]));
            
            F = f(Xphys,Yphys);
            nBasis = size(obj.Basis,3);
            U = zeros(n_cells,nBasis);
            for i=1:nBasis
                U(:,i) = (obj.massmatrix(i)).^(-1)*sum(sum(bsxfun(@times,F,permute(obj.QuadratureHR.W.*obj.BasisHR(:,:,i),[3 1 2])),3),2);
            end
        end
        
               
        function plot(obj,U,flat,logarithmic)
            if nargin<3
                flat = true;
            end
           if nargin<4
                logarithmic = false;
            end
            Uloc = sum(bsxfun(@times,permute(U,[1 3 4 2]),permute(obj.Basis,[4 1 2 3])),4);
            if logarithmic
                Uloc = log10(max(Uloc,1e-10));
            end
            
            X = reshape(obj.midpoints(obj.InteriorCells,1),[obj.ny,obj.nx]);
            Y = reshape(obj.midpoints(obj.InteriorCells,2),[obj.ny,obj.nx]);
            Uloc = reshape(Uloc,[obj.ny,obj.nx]);
            if flat
                pcolor(X,Y,Uloc);               
            else               
               surf(X,Y,Uloc)
               view(3)
           end
        end
        
         function GridData = save(obj,varargin)
            warning('off','MATLAB:structOnObject');
            GridData = struct(obj);
            GridData.Name = class(obj);
            warning('on','MATLAB:structOnObject');
            if nargin==3
                path = varargin{1}; 
                filename = varargin{2}; 
                save([path filename],'GridData')
            end
        end
        
        function load(obj,varargin)
            if nargin==2
                GridData = varargin{1}; 
            else
                path = varargin{1}; 
                filename = varargin{2}; 
                load([path filename]); %gives GridData
            end
            if strcmp(class(obj),GridData.Name)
                GridData = rmfield(GridData,'Name');
                F = fieldnames(GridData);
                for i=1:length(F)
                    obj.(F{i}) = GridData.(F{i});
                end
            else
                error(sprintf('It appears that you wanted to load the data for class "%s" into class "%s"\n',GridData.Name,class(obj))); %#ok<SPERR>
            end
            
        end
        
    end
end

