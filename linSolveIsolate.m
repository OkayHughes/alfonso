
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
    
    inter = soln.L * soln.L';
    [U, S, V] = svd(inter);
    sing = diag(S);
    sing_inv = 1./sing;
    sing_inv(sing < 1e-10) = 0;
    inver = V * diag(sing_inv) * U';
    
    Hic     = inver * c;
    HiAt    = -inver * A';
    Hirxrs  = inver * (rx+rs);
    fprintf("cond(H) = %5e\n", cond(soln.L * soln.L'))
    
    Hic_alt     = soln.L'\(soln.L\c);
    
    HiAt_alt    = -soln.L'\(soln.L\A');
    
    Hirxrs_alt  = soln.L'\(soln.L\(rx+rs));
    
    
    fprintf("Hic residual: %5d\n", norm(soln.L * soln.L' * Hic - c))
    fprintf("HiAt residual: %5d\n", norm(soln.L * soln.L' * HiAt + A'))
    fprintf("Hirxrs residual: %5d\n", norm(soln.L * soln.L' * Hirxrs - (rx + rs)))
    fprintf("Hic_alt residual: %5d\n", norm(soln.L * soln.L' * Hic_alt - c))
    fprintf("HiAt_alt residual: %5d\n", norm(soln.L * soln.L' * HiAt_alt + A'))
    fprintf("Hirxrs_alt residual: %5d\n", norm(soln.L * soln.L' * Hirxrs_alt - (rx + rs)))
    

%     f = figure('visible','off');
%     eigens = eig(soln.L * soln.L');
%     plot(sort(eigens))
%     saveas(f,sprintf('plots/hess_%d', figcount),'png')
    
    fprintf("cond(At) = %5e\n", cond(A'));
    fprintf("cond(HiAt) = %5e\n", cond(HiAt));
    fprintf("cond(Hmat) = %5e\n", cond([HiAt, Hic]));
    oof = [A; -c']*[HiAt, Hic]; 
    fprintf("cond(hmat product) = %5e\n", cond(oof))
    LHSdydtau   = [zeros(m), -b; b', soln.mu/soln.tau^2] - [A; -c']*[HiAt, Hic]/soln.mu;
%     f = figure('visible','off');
%     global figcount;
%     figcount = figcount + 1;
%     eigens = eig(LHSdydtau);
%     plot(sort(eigens))
%     saveas(f,sprintf('plots/LHS_%d', figcount),'png')
%     fprintf("cond(LHS) = %5e\n", cond(LHSdydtau));
    RHSdydtau   = [ry; rtau+rkappa] - [A; -c']*Hirxrs/soln.mu;
    dydtau      = LHSdydtau\RHSdydtau;
    fprintf("LHS residual: %5d\n\n", norm(LHSdydtau * dydtau - RHSdydtau))
    dx          = (Hirxrs - [HiAt, Hic]*dydtau)/soln.mu;

    delta               = zeros(m+2*n+2, 1);
    delta(1:m)          = dydtau(1:m);
    delta(m+n+1)        = dydtau(m+1);
    delta(m+(1:n))      = dx;
    delta(m+n+1+(1:n))  = -rx - [A', -c]*dydtau;
    delta(end)          = -rtau + b'*dydtau(1:m) - c'*dx;
    
return