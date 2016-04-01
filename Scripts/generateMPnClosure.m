function [P, Q, S] = generateMPnClosure(n)
%n =6; %Order
% syms mue phi real
% Omegax = sqrt(1-mue.^2).*cos(phi);
% Omegay = sqrt(1-mue.^2).*sin(phi);

pi = sym('pi');
M = 2*n^2+2*n+1; %Number of moments
nq = 1/2.*n.^2-1/2.*n; %Number of quarter moments per quadrant
A = sym(zeros(M,M)); %A(i,j) = intS(b(i)*b(j))
LB = A; %LB(i,j) = intS(b(j)*D_Omega b(i))
[I,J] = meshgrid(1:n-1);
ind = I+J<=n;
I = I(ind);
J = J(ind);
i = 1:n;
intSxp = sym(2*pi./(i+1));
intSxm = sym((-1).^i*2*pi./(i+1));
intSyp = sym(2*pi./(i+1));
intSym = sym((-1).^i*2*pi./(i+1));

A(1,1:4*n+1) = [1,intSxp,intSxm,intSyp,intSym];
A(1:4*n+1,1) = [1,intSxp,intSxm,intSyp,intSym]';


m = i';

A(1,1:4.*n+1) = [4.*pi,2.*pi./(i+1),(-1).^i.*2.*pi./(i+1),2.*pi./(i+1),(-1).^i.*2.*pi./(i+1)]; %E
LB(1,:) = 0;
%alpha
A(2:n+1,1) = 2.*pi./(m+1); %Fxp
A(n+2:2.*n+1,1) = (-1).^m.*2.*pi./(m+1); %Fxm
A(2.*n+2:3.*n+1,1) = 2.*pi./(m+1); %#ok<*BDSCA> %Fyp
A(3.*n+2:4.*n+1,1) = (-1).^m.*2.*pi./(m+1); %Fym

%beta
[i,m] = meshgrid(1:n);
A(2:n+1,2:4.*n+1) = [2.*pi./(i+m+1),0.*m,1/2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),(-1).^i./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2)];
A(n+2:2.*n+1,2:4.*n+1) = [0.*m,(-1).^(i+m).*2.*pi./(i+m+1),(-1).^m./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),(-1).^(m+i)./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2)];
A(2.*n+2:3.*n+1,2:4.*n+1) = [1/2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),(-1).^i./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),2.*pi./(i+m+1),0.*m];
A(3.*n+2:4.*n+1,2:4.*n+1) = [(-1).^m./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),(-1).^(m+i)./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),0.*m,(-1).^(i+m).*2.*pi./(i+m+1)];

cnt = 0;
m = 1:n;
%gamma
for i=2:n
    for k=1:i-1
        %E
        val = 1/2.*B(1/2,1+i/2).*B((k+1)/2,(i-k+1)/2);
        A(1,4.*n+2+cnt) = val;
        A(1,4.*n+2+nq+cnt) = val.*(-1).^k;
        A(1,4.*n+2+2.*nq+cnt) = val.*(-1).^i;
        A(1,4.*n+2+3.*nq+cnt) = val.*(-1).^(i-k);
        
        %Fxp
        val = 1/2.*B(1/2,1+(i+m)/2).*B((k+m+1)/2,(i-k+1)/2);
        valE = 1/2.*B(1/2,1+(i+0)/2).*B((k+0+1)/2,(i-k+1)/2);
        A(2:n+1,4.*n+2+cnt) = val';
        A(2:n+1,4.*n+2+3.*nq+cnt) = (-1).^(i-k).*val';
        
        %Fxm
        %val = 1/2.*B(1/2,1+(i+m)/2).*B((k+m+1)/2,(i-k+1)/2); %the same as
        %%above
        A(n+2:2.*n+1,4.*n+2+nq+cnt) = (-1).^(k+m').*val';
        A(n+2:2.*n+1,4.*n+2+2.*nq+cnt) = (-1).^(i+m').*val';
        
        %Fyp
        val = 1/2.*B(1/2,1+(i+m)/2).*B((k+1)/2,(i-k+m+1)/2);
        valE = 1/2.*B(1/2,1+(i)/2).*B((k+1)/2,(i-k+1)/2);
        A(2.*n+2:3.*n+1,4.*n+2+cnt) = val';
        A(2.*n+2:3.*n+1,4.*n+2+nq+cnt) = (-1).^(k).*val';
        
        %Fyp
        %val = 1/2.*B(1/2,1+(i+m)/2).*B((k+1)/2,(i-k+m+1)/2);
        A(3.*n+2:4.*n+1,4.*n+2+2.*nq+cnt) = (-1).^(i+m').*val';
        A(3.*n+2:4.*n+1,4.*n+2+3.*nq+cnt) = (-1).^(i-k+m').*val';
        
        cnt = cnt+1;
        
    end
end
cntx = 0;
%Quarter Moments
for m=2:n
    for K=1:m-1
        i = 1:n;
        A(1+4.*n+1+cntx,1:4.*n+1) = [1/2.*B(1/2,1+1/2.*m).*B((K+1)/2,(m-K+1)/2),  1/2.*B(1/2,1+1/2.*(i+m)).*B((i+K+1)/2,(m-K+1)/2),  0.*i,  1/2.*B(1/2,1+1/2.*(i+m)).*B((K+1)/2,(i+m-K+1)/2),  0.*i];
        A(1+4.*n+nq+1+cntx,1:4.*n+1) = [(-1).^K/2.*B(1/2,1+1/2.*m).*B((K+1)/2,(m-K+1)/2),  0.*i,  (-1).^(i+K)./2.*B(1/2,1+1/2.*(i+m)).*B((i+K+1)/2,(m-K+1)/2),    (-1).^K./2.*B(1/2,1+1/2.*(i+m)).*B((K+1)/2,(i+m-K+1)/2),   0.*i];
        A(1+4.*n+2.*nq+1+cntx,1:4.*n+1) = [(-1).^m/2.*B(1/2,1+1/2.*m).*B((K+1)/2,(m-K+1)/2),  0.*i,  (-1).^(i+m)./2.*B(1/2,1+1/2.*(i+m)).*B((i+K+1)/2,(m-K+1)/2),   0.*i,    (-1).^(i+m)./2.*B(1/2,1+1/2.*(i+m)).*B((K+1)/2,(i+m-K+1)/2)];
        A(1+4.*n+3.*nq+1+cntx,1:4.*n+1) = [(-1).^(m-K)/2.*B(1/2,1+1/2.*m).*B((K+1)/2,(m-K+1)/2),   (-1).^(m-K)./2.*B(1/2,1+1/2.*(i+m)).*B((i+K+1)/2,(m-K+1)/2),   0.*i,    0.*i,   (-1).^(i+m-K)./2.*B(1/2,1+1/2.*(i+m)).*B((K+1)/2,(i+m-K+1)/2)];
        
        cnty = 0;
        for i=2:n
            for k=1:i-1
                val = 1/2.*B(1/2,1+1/2.*(i+m)).*B((k+K+1)/2,(i-k+m-K+1)/2);
                A(1+4.*n+1+cntx,4.*n+2+cnty) = val;
                A(1+4.*n+nq+1+cntx,4.*n+2+nq+cnty) = (-1).^(k+K).*val;
                A(1+4.*n+2.*nq+1+cntx,4.*n+2+2.*nq+cnty) = (-1).^(i+m).*val;
                A(1+4.*n+3.*nq+1+cntx,4.*n+2+3.*nq+cnty) = (-1).^(i+m-K-k).*val;
                
                cnty = cnty+1;
            end
            
        end
        
        cntx = cntx+1;
        
    end
end
disp('Matrix done');
n2 = n+1;
nq2 = 1/2.*n2.^2-1/2.*n2;

%% Closure
RHS = sym(zeros(4.*(1+nq2-nq),M));

% i = sym(1):n;
% m = sym(n2);
i = 1:n;
m = n2;

RHS(1:4,1) = [2.*pi./(m+1);(-1).^m.*2.*pi./(m+1); 2.*pi./(m+1);(-1).^m.*2.*pi./(m+1)]; %alpha

RHS(1,2:4.*n+1) = [2.*pi./(i+m+1),0.*i,1/2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),(-1).^i./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2)];
RHS(2,2:4.*n+1) = [0.*i,(-1).^(i+m).*2.*pi./(i+m+1),(-1).^m./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),(-1).^(m+i)./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2)];
RHS(3,2:4.*n+1) = [1/2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),(-1).^i./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),2.*pi./(i+m+1),0.*i];
RHS(4,2:4.*n+1) = [(-1).^m./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),(-1).^(m+i)./2.*B(1/2,1+1/2.*(m+i)).*B((m+1)/2,(i+1)/2),0.*i,(-1).^(i+m).*2.*pi./(i+m+1)];


cnt = 0;
%gamma
for i=2:n
    for k=1:i-1
        %Fxp
        val = 1/2.*B(1/2,1+(i+m)/2).*B((k+m+1)/2,(i-k+1)/2);
        RHS(1,4.*n+2+cnt) = val';
        RHS(1,4.*n+2+3.*nq+cnt) = (-1).^(i-k).*val';
        
        %Fxm
        %val = 1/2.*B(1/2,1+(i+m)/2).*B((k+m+1)/2,(i-k+1)/2); %the same as
        %%above
        RHS(2,4.*n+2+nq+cnt) = (-1).^(k+m').*val';
        RHS(2,4.*n+2+2.*nq+cnt) = (-1).^(i+m').*val';
        
        %Fyp
        val = 1/2.*B(1/2,1+(i+m)/2).*B((k+1)/2,(i-k+m+1)/2);
        RHS(3,4.*n+2+cnt) = val';
        RHS(3,4.*n+2+nq+cnt) = (-1).^(k).*val';
        
        %Fyp
        %val = 1/2.*B(1/2,1+(i+m)/2).*B((k+1)/2,(i-k+m+1)/2);
        RHS(4,4.*n+2+2.*nq+cnt) = (-1).^(i+m').*val';
        RHS(4,4.*n+2+3.*nq+cnt) = (-1).^(i-k+m').*val';
        
        cnt = cnt+1;
        
    end
end

cntx = 0;
%Quarter Moments
for m=n2:n2
    for K=1:m-1
        i = 1:n;
        RHS(4+1+cntx,1:4.*n+1) = [1/2.*B(1/2,1+1/2.*m).*B((K+1)/2,(m-K+1)/2),  1/2.*B(1/2,1+1/2.*(i+m)).*B((i+K+1)/2,(m-K+1)/2),  0.*i,  1/2.*B(1/2,1+1/2.*(i+m)).*B((K+1)/2,(i+m-K+1)/2),  0.*i];
        RHS(4+(nq2-nq)+1+cntx,1:4.*n+1) = [(-1).^K/2.*B(1/2,1+1/2.*m).*B((K+1)/2,(m-K+1)/2),  0.*i,  (-1).^(i+K)./2.*B(1/2,1+1/2.*(i+m)).*B((i+K+1)/2,(m-K+1)/2),    (-1).^K./2.*B(1/2,1+1/2.*(i+m)).*B((K+1)/2,(i+m-K+1)/2),   0.*i];
        RHS(4+2.*(nq2-nq)+1+cntx,1:4.*n+1) = [(-1).^m/2.*B(1/2,1+1/2.*m).*B((K+1)/2,(m-K+1)/2),  0.*i,  (-1).^(i+m)./2.*B(1/2,1+1/2.*(i+m)).*B((i+K+1)/2,(m-K+1)/2),   0.*i,    (-1).^(i+m)./2.*B(1/2,1+1/2.*(i+m)).*B((K+1)/2,(i+m-K+1)/2)];
        RHS(4+3.*(nq2-nq)+1+cntx,1:4.*n+1) = [(-1).^(m-K)/2.*B(1/2,1+1/2.*m).*B((K+1)/2,(m-K+1)/2),   (-1).^(m-K)./2.*B(1/2,1+1/2.*(i+m)).*B((i+K+1)/2,(m-K+1)/2),   0.*i,    0.*i,   (-1).^(i+m-K)./2.*B(1/2,1+1/2.*(i+m)).*B((K+1)/2,(i+m-K+1)/2)];
        
        
        cnty = 0;
        for i=2:n
            for k=1:i-1
                val = 1/2.*B(1/2,1+1/2.*(i+m)).*B((k+K+1)/2,(i-k+m-K+1)/2);
                RHS(4+1+cntx,4.*n+2+cnty) = val;
                RHS(4+(nq2-nq)+1+cntx,4.*n+2+nq+cnty) = (-1).^(k+K).*val;
                RHS(4+2.*(nq2-nq)+1+cntx,4.*n+2+2.*nq+cnty) = (-1).^(i+m).*val;
                RHS(4+3.*(nq2-nq)+1+cntx,4.*n+2+3.*nq+cnty) = (-1).^(i+m-K-k).*val;
                cnty = cnty+1;
            end
            
        end
        
        cntx = cntx+1;
        
    end
end


%% Laplace Beltrami
momentsdd_pihalf = sym(zeros(length(A),1));%intS(b*delta(phi-pi/2)/sqrt(1-mu^2))
momentsdd_pihalf(1) = B(1/2,1-1/2);
momentsdd_threepihalf = momentsdd_pihalf;%intS(b*delta(phi-3*pi/2)/sqrt(1-mu^2))
momentsdd_pi = momentsdd_pihalf;%intS(b*delta(phi-pi)/sqrt(1-mu^2))
momentsdd_twopi = momentsdd_pihalf;%intS(b*delta(phi-2*pi)/sqrt(1-mu^2))
i = 1:n;
momentsdd_pihalf(1+i) = B(1/2,1+1/2*(i-1)).*cos(pi/2).^i; %intSxp(mxp*delta(...)/sqrt(..)) = 0
momentsdd_pihalf(1+n+i) = 0.5*B(1/2,1+1/2*(i-1)).*cos(pi/2).^i; %intSxp(mxm*delta(...)/sqrt(..)) = 0 %0.5 wegen dirac gewicht???
momentsdd_pihalf(1+2*n+i) = B(1/2,1+1/2*(i-1)).*sin(pi/2).^i; %intSxp(myp*delta(...)/sqrt(..))
momentsdd_pihalf(1+3*n+i) = 0*B(1/2,1+1/2*(i-1)).*sin(pi/2).^i; %intSxp(mym*delta(...)/sqrt(..)) = 0 % pi/2 nicht in Sym

%Bei den anderen entsprechend auch nur ein Eintrag ungleich 0:
momentsdd_pi(1+n+i) = B(1/2,1+1/2*(i-1)).*cos(pi).^i; %intSyp(mxm*delta(...)/sqrt(..))
momentsdd_threepihalf(1+3*n+i) = B(1/2,1+1/2*(i-1)).*sin(3*pi/2).^i; %intSxm(myp*delta(...)/sqrt(..)) = intSxp(myp*delta(...)/sqrt(..))
momentsdd_twopi(1+i) = B(1/2,1+1/2*(i-1)).*cos(2*pi).^i;
%Quarter Moments sind immer 0 da an den Grenzen entweder cos oder sin 0
%sind.

moments_1sxp = sym(zeros(length(A),1));
moments_1sxm = sym(zeros(length(A),1));
moments_1syp = sym(zeros(length(A),1));
moments_1sym = sym(zeros(length(A),1));
moments_1sxp(1:4*n+1) = [2*pi,2*pi./(i+1),0*i,1/2*B(1/2,1+1/2*i).*B(1/2,1/2+1/2*i),(-1).^i/2.*B(1/2,1+1/2*i).*B(1/2,1/2+1/2*i)];
moments_1sxp(4*n+2:end) = A(1,4*n+2:end).*[ones(1,nq),zeros(1,nq),zeros(1,nq),ones(1,nq)];
moments_1sxm(1:4*n+1) = [2*pi,0*i,(-1).^i.*2*pi./(i+1),1/2*B(1/2,1+1/2*i).*B(1/2,1/2+1/2*i),(-1).^i/2.*B(1/2,1+1/2*i).*B(1/2,1/2+1/2*i)];
moments_1sxm(4*n+2:end) = A(1,4*n+2:end).*[zeros(1,nq),ones(1,nq),ones(1,nq),zeros(1,nq)];
moments_1syp(1:4*n+1) = [2*pi,1/2*B(1/2,1+1/2*i).*B(1/2,1/2+1/2*i),(-1).^i/2.*B(1/2,1+1/2*i).*B(1/2,1/2+1/2*i),2*pi./(i+1),0*i];
moments_1syp(4*n+2:end) = A(1,4*n+2:end).*[ones(1,nq),ones(1,nq),zeros(1,nq),zeros(1,nq)];
moments_1sym(1:4*n+1) = [2*pi,1/2*B(1/2,1+1/2*i).*B(1/2,1/2+1/2*i),(-1).^i/2.*B(1/2,1+1/2*i).*B(1/2,1/2+1/2*i),0*i,(-1).^i.*2*pi./(i+1)];
moments_1sym(4*n+2:end) = A(1,4*n+2:end).*[zeros(1,nq),zeros(1,nq),ones(1,nq),ones(1,nq)];
for l=3:n
    LB(:,l+1) = -l*(l+1)*A(:,l+1)+l*(l-1)*A(:,l-1);
    LB(:,l+1+n) = -l*(l+1)*A(:,l+1+n)+l*(l-1)*A(:,l-1+n);
    LB(:,l+1+2*n) = -l*(l+1)*A(:,l+1+2*n)+l*(l-1)*A(:,l-1+2*n);
    LB(:,l+1+3*n) = -l*(l+1)*A(:,l+1+3*n)+l*(l-1)*A(:,l-1+3*n);
end
l=1;
    LB(:,l+1) = -l*(l+1)*A(:,l+1)+momentsdd_pihalf+momentsdd_threepihalf;
    LB(:,l+1+n) = -l*(l+1)*A(:,l+1+n)-momentsdd_pihalf-momentsdd_threepihalf;
    LB(:,l+1+2*n) = -l*(l+1)*A(:,l+1+2*n)+momentsdd_twopi+momentsdd_pi;
    LB(:,l+1+3*n) = -l*(l+1)*A(:,l+1+3*n)-momentsdd_twopi-momentsdd_pi;
if n>=2
    l=2;

    LB(:,l+1) = -l*(l+1)*A(:,l+1)+l*(l-1)*moments_1sxp;
    LB(:,l+1+n) = -l*(l+1)*A(:,l+1+n)+l*(l-1)*moments_1sxm;
    LB(:,l+1+2*n) = -l*(l+1)*A(:,l+1+2*n)+l*(l-1)*moments_1syp;
    LB(:,l+1+3*n) = -l*(l+1)*A(:,l+1+3*n)+l*(l-1)*moments_1sym;
end
%Quarter moments
cntx = 0;
for m=sym(2:n)
    for K=sym(1:m-1)
        cnty = 0;
                for m2=2:n
                    for K2=1:m2-1
                        val = 1/2.*B(1/2,1+1/2.*(m+m2)).*B((K+K2+1)/2,(m-K+m2-K2+1)/2); %int(m^mK*m^(ik))
                        val2 = 1/2.*B(1/2,1+1/2.*(m+m2-2)).*B((K+K2+1)/2,(m-K+m2-2-K2+1)/2); %int(m^(m,K-2)*m^(ik))=int(m^(m-2,K)*m^(ik))
                        val3 = 1/2.*B(1/2,1+1/2.*(m+m2-2)).*B((K+K2-1)/2,(m-K+m2-K2+1)/2); %int(m^(m,K-2)*m^(ik))=int(m^(m-2,K)*m^(ik))
                        %diracs sind immer 0
                        LB(1+4.*n+1+cntx,1+4.*n+1+cnty) = -(m+m^2)*val + val2*(m-K)*(m-K-1)+val3*K*(K-1);
                        LB(1+4.*n+nq+1+cntx,1+4.*n+nq+1+cnty) = (-1).^(K2+K)*(-(m+m^2)*val + val2*(m-K)*(m-K-1)+val3*K*(K-1));
                        LB(1+4.*n+2.*nq+1+cntx,1+4.*n+2.*nq+1+cnty) = (-1).^(m2+m)*(-(m+m^2)*val + val2*(m-K)*(m-K-1)+val3*K*(K-1));
                        LB(1+4.*n+3.*nq+1+cntx,1+4.*n+3.*nq+1+cnty) = (-1).^(m2+m-K-K2)*(-(m+m^2)*val + val2*(m-K)*(m-K-1)+val3*K*(K-1));
                        cnty = cnty+1;
                    end
                end
        cntx = cntx+1;
        
    end
end
LB = tril(LB)-diag(diag(LB))+tril(LB)'; %make it symmetric since not everything is computed above


%% Final
%tic, P2 = double(A\RHS'); toc
%tic, P = double(RHS.*inv(A)); toc
Ad = double(A);
LBd = double(LB);

Q = Ad(1,:)';
% P = double(A)\double(RHS');
P = linsolve(Ad,double(RHS'));
%P-P2'
P = P';
% Du = transpose(double(A)\double(D'));
%Du = transpose(linsolve(Ad,double(RHS2')));
S = mrdivide(LBd,Ad);
%Ainv = inv(A);
%Du2 = simple(RHS2*Ainv);
end

function y = B(u,v)

u(u==0) = 1e-7;
v(v==0) = 1e-7;
y = beta(u,v);
end
