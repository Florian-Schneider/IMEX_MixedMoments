function du = DifferentialOperator( obj,U,t)
%DIFFERENTIALOPERATOR Heart of the kinetic scheme. Assembles the
%explicit part of the
%discretized differential operator.
if nargin<2
    U = obj.u;
end

Q = obj.Case.getQ(obj.Grid,obj.MomentModel,t);
[BoundaryData] = obj.Case.getBC(obj.Grid,obj.MomentModel,t);
UB = cat(1,U,BoundaryData);
UB = permute(UB,[3 1 4 5 2]);
UVal = sum(bsxfun(@times,UB,permute(obj.Grid.Basis,[4 5 1 2 3])),5);

Q = permute(Q,[3 1 4 5 2]);
QVal = sum(bsxfun(@times,Q,permute(obj.Grid.Basis,[4 5 1 2 3])),5);

scatteringtrue = false; %only in the implicit part, RHS is turned off in the explicit part
[FxpB,FxmB,FypB,FymB] = obj.MomentModel.Closure(UVal,obj,scatteringtrue);
Fxp = FxpB(:,obj.Grid.InteriorCells,:,:);
Fxm = FxmB(:,obj.Grid.InteriorCells,:,:);
Fyp = FypB(:,obj.Grid.InteriorCells,:,:);
Fym = FymB(:,obj.Grid.InteriorCells,:,:);


FxpLeftN = FxpB(:,obj.Grid.LeftNeighbor,obj.Grid.RightEdge_Quadrature_Points);
FxpRight = Fxp(:,:,obj.Grid.RightEdge_Quadrature_Points);
FxmLeft = Fxm(:,:,obj.Grid.LeftEdge_Quadrature_Points);
FxmRightN = FxmB(:,obj.Grid.RightNeighbor,obj.Grid.LeftEdge_Quadrature_Points);

FypLowerN = FypB(:,obj.Grid.LowerNeighbor,obj.Grid.UpperEdge_Quadrature_Points);
FypUpper = Fyp(:,:,obj.Grid.UpperEdge_Quadrature_Points);
FymLower = Fym(:,:,obj.Grid.LowerEdge_Quadrature_Points);
FymUpperN = FymB(:,obj.Grid.UpperNeighbor,obj.Grid.LowerEdge_Quadrature_Points);

HxLeft = -FxpLeftN-FxmLeft; %Kinetic flux (Upwind)
HxRight = FxpRight+FxmRightN;
HyLower = -FypLowerN-FymLower;
HyUpper = FypUpper+FymUpperN;


W = obj.Grid.Quadrature.W;
WLeft = obj.Grid.LeftEdge_Quadrature_Weights;
WRight = obj.Grid.RightEdge_Quadrature_Weights;
WLower = obj.Grid.LowerEdge_Quadrature_Weights;
WUpper = obj.Grid.UpperEdge_Quadrature_Weights;

InterfaceIntegral = sum(bsxfun(@times,permute(HxLeft,[1 2 4 3]),permute(bsxfun(@times,obj.Grid.BasisLeft,WLeft'),[3 4 1 2])),4) ...
    +sum(bsxfun(@times,permute(HxRight,[1 2 4 3]),permute(bsxfun(@times,obj.Grid.BasisRight,WRight'),[3 4 1 2])),4) ...
    +sum(bsxfun(@times,permute(HyLower,[1 2 4 3]),permute(bsxfun(@times,obj.Grid.BasisLower,WLower'),[3 4 1 2])),4) ...
    +sum(bsxfun(@times,permute(HyUpper,[1 2 4 3]),permute(bsxfun(@times,obj.Grid.BasisUpper,WUpper'),[3 4 1 2])),4);
SourceIntegral = sum(sum(bsxfun(@times,bsxfun(@times,permute(W,[3 4 5 1 2]),permute(QVal,[1 2 5 3 4])),permute(obj.Grid.Basis,[4 5 3 1 2])),5),4);                   
du = bsxfun(@times,permute(1./obj.Grid.massmatrix,[2 3 1]),-InterfaceIntegral/obj.Grid.cell_volumes +SourceIntegral);

du = permute(du,[2 3 1]);
end


