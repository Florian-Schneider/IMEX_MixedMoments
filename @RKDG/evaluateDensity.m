function v = evaluateDensity(obj,x,y,t)

v = zeros(length(x),length(t));
u0 = permute(obj.MomentModel.Density(permute(obj.u_Frame,[3 1 2 4])),[4 2 3 1])/max(obj.MomentModel.BasisFunctionsAtQuadratureSqueezed(1,:));

DT =delaunayTriangulation(obj.Grid.midpoints(obj.Grid.InteriorCells,:));
inds = DT.nearestNeighbor([x,y]);
X = x;
Y = y;
for i=1:length(x)
    X(i) = (x(i)-obj.Grid.midpoints(inds(i),1))/obj.Grid.dx;
    Y(i) = (y(i)-obj.Grid.midpoints(inds(i),2))/obj.Grid.dy;
end

I = obj.Grid.BasisOrder(:,1);
J = obj.Grid.BasisOrder(:,2);
p = zeros([1,size(X,1),length(I)]);
for i=1:length(I)
    [px] = legPoly(2*X,I(i)); %legendre polynomials
    [py] = legPoly(2*Y,J(i)); %legendre polynomials
    p(1,:,i) = px.*py;
end
DensityEv = sum(bsxfun(@times,u0(:,inds,:),p),3);

for i=1:length(t)
    
    Density = permute(interp1(obj.t_Frame',DensityEv,t(i)),[2 3 1]);
    v(:,i) = Density;
end

end