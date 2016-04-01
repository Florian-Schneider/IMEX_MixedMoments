function [ xw ] = generateHalfRangeChebyshevPolysNum( n )
syms x real
k = 0:2*n-1;
%b = sym(1./(4*(4-k.^-2)))';
b = (1./(4*(4-k.^-2)))';
b(1) = 1;
a = 1/2+0*b;
ab = [a,b];
P{1} = sym(1);
P{2} = (x-a(1))*P{1}-b(1)*0;
m = (int(sqrt(1-x^2)^-1*x.^(0:2*n-1),x,0,1));
m2 = m;
i=1;
% p = sym2poly(P{i});
% m2(i) = fliplr(p)*m(1:i)';
p0 = coeffs(P{i});
m2(i) = p0*m(1:i)';
i=2;
% p = sym2poly(P{i});
% m2(i) = fliplr(p)*m(1:i)';
p1 = coeffs(P{i});
m2(i) = p1*m(1:i)';
for i=3:2*n
    p2 = [-a(i)*p1(1)-b(i-1)*p0(1),p1(1:end-2)-a(i)*p1(2:end-1)-b(i-1)*p0(2:end),p1(end-1)-a(i)*p1(end),p1(end)];
    p0 = p1;
    p1 = p2;
    m2(i) = p2*m(1:i)';
end

[ab2]= chebyshev(n,m2,ab);
% a = ab2(:,1); b = ab2(:,2);
% P2{1} = 1;
% P2{2} = (x-a(1))*P2{1}-b(1)*0;
% m(1) = int(P2{1}/sqrt(1-x^2),x,0,1);
% m(2) = int(P2{2}/sqrt(1-x^2),x,0,1);
% for i=3:n+1
%     P2{i} = (x-a(i-1))*P2{i-1}-b(i-1)*P2{i-2};
%     m(i) = int(P2{i}/sqrt(1-x^2),x,0,1);
% end

xw = double(gauss(n,ab2));
end

