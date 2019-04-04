
function [delta, probData] = linSolve3(soln, probData, RHS)
% This method implements a block solver for the Newton system.
% --------------------------------------------------------------------------
% USAGE of "linSolve3"
% delta = linSolve3(soln, probData, RHS)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% RHS:          right-hand side of the Newton system
%
% OUTPUT
% delta:	computed Newton direction
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------
    
    A = probData.A;
    b = probData.b;
    c = probData.c;
    [m, n] = size(A);
    
    ry     = RHS(1:m);
    rx     = RHS(m+(1:n));
    rtau   = RHS(m+n+1);
    rs     = RHS(m+n+1+(1:n));
    rkappa = RHS(end);
    
    Hic     = soln.L'\(soln.L\c);
    HiAt    = -soln.L'\(soln.L\A');
    Hirxrs  = soln.L'\(soln.L\(rx+rs));
    fprintf("cond(H) = %5e", cond(soln.L * soln.L'))
    fprintf("cond(At) = %5e\n", cond(A'));
    fprintf("cond(HiAt) = %5e\n", cond(HiAt));
    fprintf("cond(Hmat) = %5e\n", cond([HiAt, Hic]));
    oof = [A; -c']*[HiAt, Hic]; 
    fprintf("cond(hmat product) = %5e\n", cond(oof))
    LHSdydtau   = [zeros(m), -b; b', soln.mu/soln.tau^2] - [A; -c']*[HiAt, Hic]/soln.mu;
    fprintf("cond(LHS) = %5e\n", cond(LHSdydtau));
    RHSdydtau   = [ry; rtau+rkappa] - [A; -c']*Hirxrs/soln.mu;
    dydtau      = LHSdydtau\RHSdydtau;
    dx          = (Hirxrs - [HiAt, Hic]*dydtau)/soln.mu;

    delta               = zeros(m+2*n+2, 1);
    delta(1:m)          = dydtau(1:m);
    delta(m+n+1)        = dydtau(m+1);
    delta(m+(1:n))      = dx;
    delta(m+n+1+(1:n))  = -rx - [A', -c]*dydtau;
    delta(end)          = -rtau + b'*dydtau(1:m) - c'*dx;
    
return