function [ Mx,My ] = Flux( N)
%PN Summary of this function goes here
%   Detailed explanation goes here
cnt = 1;
nmom = N^2+2*N+1;
I = zeros(N^2+2*N+1,2);
Ieven = zeros(N^2+2*N+1,1);
for l=0:N
    for m=-l:l
        I(cnt,:) = [l,m]; %Sorting of complex spherical harmonics
        if mod(l+m,2)==0
            Ieven(cnt)=1;
        end
        cnt = cnt+1;
    end
end

Mxc = zeros(nmom,nmom);
Myc = Mxc;

for i=1:nmom
    l = I(i,1); m = I(i,2);
    [~,j1] = ismember([l-1,m-1],I,'rows');
    if j1>0
        Mxc(i,j1) = -c(l-1,m-1);
        Myc(i,j1) = c(l-1,m-1);
    end
    [~,j2] = ismember([l+1,m-1],I,'rows');
    if j2>0
        Mxc(i,j2) = d(l+1,m-1);
        Myc(i,j2) = -d(l+1,m-1);
    end
    [~,j3] = ismember([l-1,m+1],I,'rows');
    if j3>0
        Mxc(i,j3) = e(l-1,m+1);
        Myc(i,j3) = e(l-1,m+1);
    end
    [~,j4] = ismember([l+1,m+1],I,'rows');
    if j4>0
        Mxc(i,j4) = -f(l+1,m+1);
        Myc(i,j4) = -f(l+1,m+1);
    end
end
Mxc = 1/2*Mxc;
Myc = 1i/2*Myc;
% for i=nmom:-1:1
%      l = I(i,1); m = I(i,2);
%      if mod(l+m,2)==1
%          Mxc(i,:) = [];
%          Mxc(:,i) = [];
%          Myc(i,:) = [];
%          Myc(:,i) = [];
%      end
% end

for l=0:N
    C{l+1} = zeros(2*l+1);
    shift = l+1;
    for m = -l:l
        C{l+1}(abs(m)+shift,abs(m)+shift) = (-1)^abs(m)/sqrt(2);
        C{l+1}(abs(m)+shift,-abs(m)+shift) = 1/sqrt(2);
        C{l+1}(-abs(m)+shift,abs(m)+shift) = -1i*(-1)^abs(m)/sqrt(2);
        C{l+1}(-abs(m)+shift,-abs(m)+shift) = 1i/sqrt(2);
        C{l+1}(0+shift,0+shift) = 1;
    end
%     for m=l:-1:-l
%         if mod(l+m,2)==1
%             C{l+1}(m+shift,:) = [];
%             C{l+1}(:,m+shift) = [];
%         end
%     end
end

S = blkdiag(C{:});
Mx = S*Mxc*S';
My = -S*Myc*S'; %don't know where this - comes from, but then it is correct
Mx(logical(1-Ieven),:) = [];
My(logical(1-Ieven),:) = [];
Mx(:,logical(1-Ieven)) = [];
My(:,logical(1-Ieven)) = [];
% for i=nmom:-1:1
%      l = I(i,1); m = I(i,2);
%      if mod(l+m,2)==1
%          Mx(i,:) = [];
%          Mx(:,i) = [];
%          My(i,:) = [];
%          My(:,i) = [];
%      end
% end
end

function y = c(l,m)
y = sqrt((l+m+1)*(l+m+2)/(2*l+3)/(2*l+1));
end
function y = d(l,m)
y = sqrt((l-m)*(l-m-1)/(2*l+1)/(2*l-1));
end
function y = e(l,m)
y = sqrt((l-m+1)*(l-m+2)/(2*l+3)/(2*l+1));
end
function y = f(l,m)
y = sqrt((l+m)*(l+m-1)/(2*l+1)/(2*l-1));
end
