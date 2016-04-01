function [alpha, beta, T, outs] = dualAdaptivePoly(u, beta0, T0, ...
    parameters)
% dualAdaptivePoly: Wrapper for the mex-File of the dual adaptive basis
% optimization solver for the moment problem <b*exp(alpha*b)> = u
% with initial guess beta0 and basis transformation matrix T0
%
% (c) Florian Schneider (schneider@mathematik.uni-kl.de) 2016
% See licence.txt for more information about usage and redistribution

if u(1)>0
    [alpha,beta,T,outs.status,outs.iter,outs.normg,outs.r,outs.rIdx,outs.maxHcond,outs.wG] = dualAdaptivePoly_mex(u,beta0,T0,parameters.itermax,parameters.rList,parameters.rList2,parameters.phiIso,parameters.p,parameters.w,parameters.k0,parameters.tol,1);
    outs.normg = outs.normg(1:min(outs.iter,length(outs.normg)));
    outs.normg0 = outs.normg(1);
else
    outs.status = 0;
    outs.maxHcond = inf;
    alpha = 0*beta0;
    beta = beta0;
    T = T0;
    
  
end





end

