
function [solnAlpha, corrStatus] = corr(soln, probData, gH, gH_Params, myLinSolve, algParams, opts)
% This method performs a single corrector step in the algorithm.
% --------------------------------------------------------------------------
% USAGE of "corr"
% [solnAlpha, corrStatus] = corr(soln, probData, gH, gH_Params, myLinSolve, algParams, opts)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% gH:           method for computing the gradient and Hessian 
%               of the barrier function
% gH_Params:    parameters associated with the method gH
% myLinSolve:   solution method for the Newton system
% algParams:    algorithmic parameters
% opts:         algorithmic options
%
% OUTPUT
% solnAlpha:    new iterate
% corrStatus:   0 if corrector phase was not successful. 1 if corrector 
%               phase was successful.
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------
    if isfield(opts, "verbose")
        if opts.verbose
            fprintf("Entering Corrector\n")
            fprintf("==================\n\n")
        end
    end
    [m, n] = size(probData.A);
    
    RHS = zeros(m+2*n+2,1);
    RHS(m+n+2:end)  = -soln.psi;
    
    dsoln = linSolveMain(soln, probData, RHS, myLinSolve, algParams, opts);
    alpha = algParams.alphaCorr;
    corrStatus = 1;
    
    % does line search to make sure next primal iterate remains inside the cone
    for nSteps = 1:opts.maxCorrLSIters
        
        solnAlpha.y     = soln.y     + alpha*dsoln.y;
        solnAlpha.x     = soln.x     + alpha*dsoln.x;
        solnAlpha.tau   = soln.tau   + alpha*dsoln.tau;
        solnAlpha.s     = soln.s     + alpha*dsoln.s;
        solnAlpha.kappa = soln.kappa + alpha*dsoln.kappa;
        
        [solnAlpha.in, solnAlpha.g, solnAlpha.H, solnAlpha.L] = gH(solnAlpha.x, gH_Params);
        
        % terminates line search if primal iterate is inside the cone
        if solnAlpha.in
            solnAlpha.mu   = (solnAlpha.x'*solnAlpha.s + solnAlpha.tau*solnAlpha.kappa)/gH_Params.bnu;
            solnAlpha.psi  = [solnAlpha.s;solnAlpha.kappa] + solnAlpha.mu*[solnAlpha.g;-1/solnAlpha.tau];             
            break;
        end
        
        alpha = algParams.corrLSMulti*alpha;        
    end
    
    if solnAlpha.in == 0
        corrStatus  = 0; % corrector has failed
        solnAlpha   = soln;
    end
    
return