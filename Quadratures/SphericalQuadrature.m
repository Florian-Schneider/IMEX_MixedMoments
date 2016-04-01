function [ Q ] = SphericalQuadrature( order )
%SPHERICALQUADRATURE Summary of this function goes here
%   Detailed explanation goes here

% nqs = order+1+double(mod(order,2)==0);
nqs = order+1;
[mu,wmu] = lgwt(nqs,-1,1);
I = mu>=0;
mu = mu(I); %we only have even functions
wmu = wmu(I);
wmu(mu>0) = 2*wmu(mu>0);
phi = pi/nqs*((1:2*nqs) -0.5);
wphi = pi/nqs*ones(size(phi));
[mu,phi] = meshgrid(mu,phi);
[wmu,wphi] = meshgrid(wmu,wphi);
w = wmu(:).*wphi(:);
Q = [mu(:),phi(:),w]';

% BF = generateBasisFunctions(Q,order);
% int = sum(bsxfun(@times,w',BF),2);
end

function [ BasisFunctionsAtQuadratureSqueezed] = generateBasisFunctions(Q,lmax)

%
mu = Q(1,:);
phi = Q(2,:);

P = cell(lmax+1);

for l=0:lmax
    L = legendre(l,mu);
    for m=1:l+1
        P{m,l+1} = L(m,:);
    end
end
C = cell(lmax,1);
S = C;
[C{:}] = deal(cos(phi));
[S{:}] = deal(sin(phi));

for m=2:lmax
    C{m} = C{1}.*C{m-1}-S{1}.*S{m-1};
    S{m} = S{1}.*C{m-1}+C{1}.*S{m-1};
end

NumberMoments = 1/2*(lmax^2+3*lmax)+1;

BasisFunctionsAtQuadrature = zeros(NumberMoments,1,size(mu,2));
[l,m] = meshgrid(0:lmax);
N = sqrt((2*l+1)/4/pi.*factorial(max(l-m,0))./factorial(l+m));
N(m>l) = 0;
%             N(1,1) = 1;
cnt = 1;
for l=0:lmax
    for m=-l:l
        if mod(l+m,2)==0
            if m<0
                BasisFunctionsAtQuadrature(cnt,1,:) = sqrt(2)*S{-m}.*N(-m+1,l+1).*P{-m+1,l+1};
            end
            if m==0
                BasisFunctionsAtQuadrature(cnt,1,:) = N(m+1,l+1).*P{m+1,l+1};
            end
            if m>0
                BasisFunctionsAtQuadrature(cnt,1,:) = sqrt(2)*C{m}.*N(m+1,l+1).*P{m+1,l+1};
            end
            cnt = cnt+1;
        end
    end
end
BasisFunctionsAtQuadratureSqueezed = squeeze(BasisFunctionsAtQuadrature);
end

