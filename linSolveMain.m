
function dsoln = linSolveMain(soln, probData, RHS, myLinSolve, algParams, opts)
% This method sets up the Newton system and computes its solution.
% --------------------------------------------------------------------------
% USAGE of "linSolveMain"
% dsoln = linSolveMain(soln, probData, RHS, myLinSolve, algParams, opts)
% --------------------------------------------------------------------------
% INPUT
% soln:         current iterate
% probData:     data for the conic optimization problem
% RHS:          right-hand side of the Newton system
% myLinSolve:   solution method for the Newton system
% algParams:    algorithmic parameters
% opts:         algorithmic options
%
% OUTPUT
% dsoln:	computed Newton direction
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% --------------------------------------------------------------------------
    
    [m, n] = size(probData.A);
    delta  = myLinSolve(soln, probData, RHS);
    
    if opts.maxItRefineSteps > 0 
               
        % checks to see if we need to refine the solution
        if rcond(full(soln.H)) < eps            
            
            LHS                         = probData.LHS;    
            LHS(m+n+1+(1:n),m+(1:n))    = soln.mu*soln.H;
            LHS(end,m+n+1)              = soln.mu/soln.tau^2;
            
            exitFlag    = 0;
            res         = residual3p(LHS, delta, RHS);
            resNorm     = norm(res);  
            for iter = 1:opts.maxItRefineSteps
                if exitFlag; break; end;
                d           = myLinSolve(soln, probData, res);
                deltaNew	= delta - d;
                resNew      = residual3p(LHS, deltaNew, RHS);
                resNewNorm	= norm(resNew);
                    
                % stops iterative refinement if there is not enough progress
                if resNewNorm > algParams.itRefineThreshold*resNorm
                    exitFlag = 1;
                end
                % updates solution if residual norm is smaller
                if resNewNorm < resNorm
                    delta       = deltaNew;
                    res         = resNew;
                    resNorm     = resNewNorm;
                end
            end
        end
        
    end
    
    dsoln.y     = delta(1:m);
    dsoln.x     = delta(m+(1:n));
    dsoln.tau   = delta(m+n+1);
    dsoln.s     = delta(m+n+1+(1:n));
    dsoln.kappa = delta(end);
    
return
