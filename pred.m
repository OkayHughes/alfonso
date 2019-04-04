
function [solnAlpha, alpha, betaAlpha, algParams, predStatus] = pred(soln, probData, gH, gH_Params, myLinSolve, algParams, opts)
% This method performs a predictor step in the algorithm.
% --------------------------------------------------------------------------
% USAGE of "pred"
% [solnAlpha, alpha, betaAlpha, algParams, predStatus] = pred(soln, probData, gH, gH_Params, myLinSolve, algParams, opts)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% gH:           method for computing the gradient and Hessian of the barrier function
% gH_Params:	parameters associated with the method gH
% myLinSolve:   solution method for the Newton system
% algParams:    algorithmic parameters
% opts:         algorithmic options
%
% OUTPUT
% solnAlpha:    new iterate
% alpha:        predictor step size
% betaAlpha:    neighborhood parameter at the end of the predictor phase
% algParams:    algorithmic parameters
% predStatus:   0 if predictor phase was not successful. 1 if predictor 
%               phase was successful.
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------

    if isfield(opts, "verbose")
        if opts.verbose
            fprintf("Entering Corrector\n\n")
        end
    end
    
    x       = soln.x;
    tau     = soln.tau;
    s       = soln.s;
    kappa   = soln.kappa;
    y       = soln.y;
    
    A = probData.A;
    b = probData.b;
    c = probData.c;
    [m, n] = size(A);
    
    RHS = zeros(m+2*n+2,1);
    RHS(1:m+n+1)        = -[A*x - b*tau; -A'*y + c*tau - s; b'*y - c'*x - kappa];
    RHS(m+n+2:end)      = -[s; kappa];
    
    % computes Newton direction
    dsoln = linSolveMain(soln, probData, RHS, myLinSolve, algParams, opts);
    
    % line search is done even when predLineSearch == 0 to make sure next
    % iterate stays within the beta-neighborhood.
    betaAlpha       = Inf; 
    solnAlphaInNhd  = 0;
    alphaPrevOK     = -1;
    predStatus      = 1;
    
    nSteps = 0;
    alpha = algParams.alphaPred;
    
    while true
        nSteps = nSteps + 1;
        
        solnAlpha.y         = y     + alpha*dsoln.y;
        solnAlpha.x         = x     + alpha*dsoln.x;
        solnAlpha.tau       = tau   + alpha*dsoln.tau;
        solnAlpha.s         = s     + alpha*dsoln.s;
        solnAlpha.kappa     = kappa + alpha*dsoln.kappa;            

        [solnAlpha.in, solnAlpha.g, solnAlpha.H, solnAlpha.L] = gH(solnAlpha.x, gH_Params);

        % primal iterate is inside the cone
        if solnAlpha.in
            solnAlpha.mu	= (solnAlpha.x'*solnAlpha.s +...
                solnAlpha.tau*solnAlpha.kappa)/gH_Params.bnu;
            solnAlpha.psi	= [solnAlpha.s; solnAlpha.kappa] +...
                solnAlpha.mu*[solnAlpha.g; -1/solnAlpha.tau];                
            betaAlpha       = sqrt(sum((solnAlpha.L\solnAlpha.psi(1:end-1)).^2) + ...
                (solnAlpha.tau*solnAlpha.psi(end))^2)/solnAlpha.mu;
            solnAlphaInNhd  = (betaAlpha < algParams.beta);
        end
        
        % iterate is inside the beta-neighborhood
        if solnAlpha.in && solnAlphaInNhd 
            % either the previous iterate was outside the beta-neighborhood
            % or increasing alpha again will make it > 1 
            if alphaPrevOK == 0 || alpha > algParams.predLSMulti
                if opts.predLineSearch == 1; algParams.alphaPred = alpha; end;
                break;
            end
            alphaPrevOK     = 1;
            alphaPrev       = alpha;
            betaAlphaPrev   = betaAlpha;
            solnAlphaPrev   = solnAlpha;
            alpha           = alpha/algParams.predLSMulti;
        else    % iterate is outside the beta-neighborhood
            % previous iterate was in the beta-neighborhood
            if alphaPrevOK == 1
                alpha       = alphaPrev;
                betaAlpha   = betaAlphaPrev;
                solnAlpha   = solnAlphaPrev;
                if opts.predLineSearch == 1; algParams.alphaPred = alpha; end;
                break;
            end
            % last two iterates were outside the beta-neighborhood and
            % alpha is very small
            if alpha < algParams.alphaPredThreshold
                predStatus  = 0; % predictor has failed
                alpha       = 0;
                betaAlpha   = Inf; % overwritten later in alfonso
                solnAlpha   = soln;
                if opts.predLineSearch == 1; algParams.alphaPred = alpha; end;
                break;
            end
            % alphaPrev, betaAlphaPrev, solnAlphaPrev will not be used
            % given alphaPrevOK == 0
            alphaPrevOK     = 0;
            alphaPrev       = alpha;
            betaAlphaPrev   = betaAlpha;
            solnAlphaPrev   = solnAlpha;
            alpha = algParams.predLSMulti*alpha;
        end
    end
        
return