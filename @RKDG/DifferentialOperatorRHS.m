function Unew = DifferentialOperatorRHS( obj,U,t)
%DIFFERENTIALOPERATORRHS Heart of the kinetic scheme. Assembles the
%implicit part of the
%discretized differential operator.
if nargin<2
    U = obj.u;
end

sigma_a = obj.Case.getSigma_a(obj.Grid,t)';
sigma_s = obj.Case.getSigma_s(obj.Grid,t)';

scatteringtrue = norm(sigma_s(:))>0;

UB = permute(U,[3 1 2]);



Unew = UB;
[~,~,~,~,UC,JC] = obj.MomentModel.Closure(Unew,obj,scatteringtrue);
feq = Unew-UB+obj.dt*bsxfun(@times,sigma_a,Unew)-obj.dt*bsxfun(@times,sigma_s/2,UC);
converged = sqrt(sum(feq.^2,1))<=1e-7;

if ~isempty(JC) %Newton
    
    Ja = repmat(diag(ones(obj.MomentModel.NumberMoments,1)),[1 1 size(UB,2)]);
    
    
    J_s = bsxfun(@times,1+obj.dt*permute(sigma_a,[1 3 2]),Ja);
    
    
    
    
    nsteps = 0;
    nmax = 50;
    while ~all(converged)
        J = J_s-bsxfun(@times,obj.dt*permute(sigma_s/2,[1 3 2]),JC);
        du = 0*UC;
        
        for i=1:size(du,2)
            if ~converged(i)
                du(:,i) = -J(:,:,i)\feq(:,i);
            end
        end
        Unew = Unew + du;
        [~,~,~,~,UC,JC] = obj.MomentModel.Closure(Unew,obj,scatteringtrue);
        feq = Unew-UB+obj.dt*bsxfun(@times,sigma_a,Unew)-obj.dt*bsxfun(@times,sigma_s/2,UC);
        converged = sqrt(sum(feq.^2,1))<=1e-7;
        nsteps = nsteps+1;
        if nsteps > nmax
            error('Newton did not converge')
        end
    end
else
    options = optimset('fsolve');
    options.Display = 'none';
    
    for i=1:size(UB,2)
        if ~converged(i)
            %Unew(:,i) = fsolve(@(U) fun(obj,U,UB(:,i),sigma_a(i),sigma_s(i)),UB(:,i),options);        
            Unew(:,i) = UB(1,i)*brsola(UB(:,i)/UB(1,i),@(U) fun(obj,U,UB(:,i)/UB(1,i),sigma_a(i),sigma_s(i)),[1e-5,1e-5]);
        end
    end
    
    
    
end


Unew = permute(Unew,[2 3 1]);
end

function y = fun(obj,Unew,UB,sigma_a,sigma_s)
scatteringtrue = true; %Otherwise jacobian would not be empty
[~,~,~,~,UC] = obj.MomentModel.Closure(Unew,obj,scatteringtrue);
y = Unew-UB+obj.dt*bsxfun(@times,sigma_a,Unew)-obj.dt*bsxfun(@times,sigma_s/2,UC);
end

