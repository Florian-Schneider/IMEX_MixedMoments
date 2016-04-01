function [ Qpp,Qmp,Qmm,Qpm ] = SphericalQuadratureQM( order )
%SPHERICALQUADRATURE Summary of this function goes here
%   Detailed explanation goes here

nqs = (order+1+double(mod(order,2)==0))/2+1;
% load HalfRangeChebyshevPolys
% nqs = min(nqs,length(xwd));
% x = xwd{nqs}(:,1);
% w = xwd{nqs}(:,2);
xw = generateHalfRangeChebyshevPolysNum(nqs);
x = xw(:,1); w = xw(:,2);
mu = sqrt(1-x.^2);
wmu = 2*x.*w;

tw = trigauss(order,0,pi/2);
phi = tw(:,1); wphi = tw(:,2);
[phi,mu] = meshgrid(phi,mu);
[wphi,wmu] = meshgrid(wphi,wmu);
w = wmu(:).*wphi(:);
Qpp = [mu(:),phi(:),w]';
Qmp = [mu(:),phi(:)+pi/2,w]';
Qmm = [mu(:),phi(:)+pi,w]';
Qpm = [mu(:),phi(:)+3/2*pi,w]';

end