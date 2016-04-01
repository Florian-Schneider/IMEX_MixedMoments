function [ quadrature ] = LebedevQuadrature( order )
%LEBEDEVQUADRATURE Summary of this function goes here
%   Detailed explanation goes here
alloweddegree =[ 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302,...
    350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074,...
    3470, 3890, 4334, 4802, 5294, 5810 ];
allowedorder =[3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,35,41,47,53,59,65,71,77,...
    83,89,95,101,107,113,119,125,131];

[~,I] = min(abs(allowedorder-order));
if allowedorder(I)<order
    order = allowedorder(I+1);
    degree = alloweddegree(I+1);
else
    order = allowedorder(I);
    degree = alloweddegree(I);
end

quadrature = load(sprintf('lebedev_%03d.txt',order));
theta = quadrature(:,2)/360*2*pi;
phi = quadrature(:,1)/360*2*pi;
mu = cos(theta);
w = 4*pi*quadrature(:,3)';
I = mu>=0;
phi = phi(I);
mu = mu(I);
w = w(I);
w(mu>=1e-15) = 2*w(mu>=1e-15);
QP = [mu';phi'];

quadrature = [QP;w];
end

