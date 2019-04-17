% This code is an implementation of the sum-of-squares optimization approach 
% based on non-symmetric conic optimization and polynomial interpolants 
% presented in:
%
% D. Papp and S. Yildiz. Sum-of-squares optimization without semidefinite 
% programming. Available at https://arxiv.org/abs/1712.01792.
%
% The implementation formulates and solves the polynomial optimization
% problems described in the same reference. 
% -------------------------------------------------------------------------
% Copyright (C) 2018 David Papp and Sercan Yildiz.
%
% Authors:  
%          David Papp       <dpapp@ncsu.edu>
%          Sercan Yildiz    <syildiz@email.unc.edu>  
%
% Date: 06/14/2018
%
% This code has been developed and tested with Matlab R2016b.
% -------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FILE
% None.
% -------------------------------------------------------------------------

function sol = polyOpt(degree)
% This is the main method for the sum-of-squares optimization approach to 
% to the polynomial optimization problem.
% --------------------------------------------------------------------------
% USAGE of "polyOpt"
% results = polyOpt(intParams, polyName, tol)
% --------------------------------------------------------------------------
% INPUT
% intParams:        data for the interpolant basis (as generated in
%                   ChebInterval.m, for instance)
% - intParams.n:    number of arguments to the polynomials
% - intParams.d:    largest degree of polynomials to be squared
% - intParams.L:    dimension of the space of (intParams.n)-variate
%                   degree-(intParams.d) polynomials
% - intParams.U:    dimension of the space of (intParams.n)-variate
%                   degree-(2*intParams.d) polynomials
% - intParams.pts:  interpolation points. (intParams.U x intParams.n) array. 
% - intParams.w:    quadrature weights for the interpolation points
% - intParams.P0:   evaluations of a basis for the space of
%                   (intParams.n)-variate degree-(intParams.d) 
%                   polynomials at the interpolation points.
%                   (intParams.U x intParams.L) array. columns are indexed
%                   with the basis polynomials. it is assumed that the first 
%                   nchoosek(intParams.n+k,intParams.n) columns are a
%                   basis for the space of (intParams.n)-variate degree-k
%                   polynomials for k = 1,...,intParams.d.
% - intParams.P:    similar to intParams.P0
% polyName:         name of polynomial to be minimized. see the function
%                   setPolyParams below for options.
% tol:              tolerance parameter for optimization
%
% OUTPUT
% results:                  final solution and iteration statistics
% - results.nIterations:	total number of iterations
% - results.alphaPred:      predictor step size at each iteration
% - results.betaPred:       neighborhood parameter at the end of the
%                           predictor phase at each iteration
% - results.etaCorr:        neighborhood parameter at the end of the
%                           corrector phase at each iteration
% - results.mu:             complementarity gap at each iteration
% - results.x:              final value of the primal variables
% - results.s:              final value of the dual slack variables
% - results.y:              final value of the dual free variables
% - results.tau:            final value of the tau-variable
% - results.kappa:          final value of the kappa-variable
% - results.pObj:           final primal objective value
% - results.dObj:           final dual objective value
% - results.dGap:           final duality gap
% - results.cGap:           final complementarity gap
% - results.rel_dGap:       final relative duality gap   
% - results.rel_cGap:       final relative complementarity gap
% - results.pRes:           final primal residuals
% - results.dRes:           final dual residuals
% - results.pIn:            final primal infeasibility
% - results.dIn:            final dual infeasibility
% - results.rel_pIn:        final relative primal infeasibility
% - results.rel_dIn:        final relative dual infeasibility
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------
    tol = 1e-5;
    intParams = FeketeCube(1, degree);
    
    n   = intParams.n;
    d   = intParams.d;
    U   = intParams.U;
    L   = intParams.L;
    P0  = intParams.P0;
    pts = intParams.pts;
    lb = -1;
    ub = 1;
    
    % transforms points to fit the domain
    scale   = (ub-lb)/2;
    shift   = (lb+ub)/2;
    pts     = bsxfun(@plus,bsxfun(@times,pts,scale'),shift');    
    wtVals  = pts; %This is the line where the polynomial weights g_i are evaluated on the interpolation points
                   %eg: this is the weight g(x) = x
                   %This defines the set x \geq 0
    %wtVals = pts - 1; %g(x) = x - 1, x \geq -1
    
    
    % ORDER OF VARIABLES 
    % x \in WSOS_(n,2*d)^*
   
    % WEIGHTS
    % weight #0: 1
    % degree of associated SOS multiplier: 2*d
    % weight #j for j = 1,...,n: (lb_j-t_j)(ub_j-t_j)
    % degree of associated SOS multiplier: 2*d-2

    LWts = nchoosek(n+d-1,n)*ones(n,1);
    
    % PARAMETERS ASSOCIATED WITH THE METHOD gH_polyOpt
    % see input description for gH_polyOpt for details

    gH_Params.n = n;
    gH_Params.d = d;
    gH_Params.U = U;
    
    gH_Params.L     = L;
    gH_Params.LWts  = LWts;
    nu              = L + sum(LWts);
    gH_Params.bnu	= nu+1;

    % P0 has columns of Chebyshev polynomial evaluations:

  
    gH_Params.P = P0;
    % associated positive semidefinite cone constraint:
    % P0'*diag(x)*P0 >= 0
    PWts = cell(n,1);
    for j = 1:n
        PWts{j} = diag(sqrt(wtVals(:,j)))*P0(:,1:LWts(j));
        % associated positive semidefinite cone constraint: 
        % PWts{j}'*diag(x)*PWts{j} >= 0
    end
    gH_Params.PWts = PWts;
    
    % x \in WSOS_(n,2*d)^* <=> 
    % P'*diag(x)*P >= 0, PWts{1}'*diag(x)*PWts{1} >= 0, PWts{2}'*diag(x)*PWts{2} >= 0,...

    % DATA FOR THE CONIC OPTIMIZATION PROBLEM
    probData.A = ones(1,U);
    probData.b = 1;
    probData.c = pts; %This pts should be a vector of the values of the objective function on the interpolant points
                      % In this case the objective function is f(x) = x

    tic
    
    % INITIAL PRIMAL ITERATE
    x0 = ones(U,1);
    [~, g0, ~, ~] = gH_polyOpt(x0, gH_Params);
    % scaling factor for the primal problem
    rP = max((1+abs(probData.b))./(1+abs(probData.A*x0)));
    % scaling factor for the dual problem
    rD = max((1+abs(g0))./(1+abs(probData.c)));
    % initial primal iterate
    x0 = repmat(sqrt(rP*rD), U, 1);
    
    % CUSTOM ALGORITHMIC OPTIONS
    opts.optimTol = tol;
    
    % CALL TO alfonso
    sol = alfonso(probData,x0,@gH_polyOpt,gH_Params,opts);
    fprintf('alfonso is done.\n');
    
    toc
    
    % prints error statistics
    fprintf('\n');
    fprintf('FINAL:\n');
    fprintf('Relative primal infeasibility: %d\n', sol.rel_pIn);
    fprintf('Relative dual infeasibility: %d\n', sol.rel_dIn);
    fprintf('Relative duality gap: %d\n', sol.rel_dGap);
    fprintf('Relative complementarity gap: %d\n\n', sol.rel_cGap);
    
return

function [in, g, H, L] = gH_polyOpt(x, params)
% This method computes the gradient and Hessian of the barrier function for
% the polynomial optimization problem.
% --------------------------------------------------------------------------
% USAGE of "gH_polyOpt"
% [in, g, H, L] = gH_polyOpt(x, params)
% --------------------------------------------------------------------------
% INPUT
% x:                    primal iterate
% params:               parameters associated with the method gH
% - params.n:           number of arguments to the polynomials
% - params.d:           largest degree of polynomials to be squared
% - params.U:           dimension of the space of (params.n)-variate 
%                       degree-(2*params.d) polynomials
% - params.L:           dimension of the space of (params.n)-variate 
%                       degree-(params.d) polynomials.
% - params.LWts:        dimensions of the "weighted" polynomial spaces. 
%                       params.LWts(j) is the dimension of the space of
%                       (params.n)-variate degree-(params.d-1) polynomials
%                       for j = 1,...,n.
% - params.bnu:         complexity parameter of the augmented barrier (nu-bar)
% - params.P:           evaluations of the basis for the space of 
%                       (params.n)-variate degree-(params.d) polynomials
%                       at the interpolation points
% - params.PWts:        evaluations of "weighted" basis polynomials at the
%                       interpolation points. params.PWts{j} is
%                       the evaluations of the basis for the "weighted"
%                       space of (params.n)-variate degree-(params.d-1)
%                       polynomials at the interpolation points for 
%                       j = 1,...,n. the weight corresponding to 
%                       params.PWts{j} is sqrt((lb_j+t_j)(ub_j-t_j)).
%
% OUTPUT
% in:	0 if x is not in the interior of the cone. 1 if x is in the 
%       interior of the cone.
% g:    gradient of the barrier function at x
% H:    Hessian of the barrier function at x
% L:    Cholesky factor of H
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    n       = params.n;
    P       = params.P;
    PWts    = params.PWts;
    
    % ORDER OF VARIABLES 
    % x \in WSOS_(n,2*d)^*
    
    % for the weight 1
    [in, g, H] = gH_SOSWt(x,P);
    
    if in == 1
        for j = 1:n
            % for the weight (lb_j+t_j)(ub_j-t_j)
            [inWt, gWt, HWt] = gH_SOSWt(x,PWts{j});
            in  = in & inWt;
            if in == 1
                g   = g+gWt;
                H   = H+HWt;
            else
                g   = NaN;
                H   = NaN;
                break;
            end
        end
    end
    
    if in == 1        
        % checks positive semidefiniteness of H one more time.
        % H may fail Cholesky factorization even if in == 1
        % due to numerical errors in summing up HWt's above.
        [L, err] = chol(H,'lower');
        in = in & (err == 0);
    else
        L = NaN;
    end
    
return

function [in, g, H] = gH_SOSWt(x, P)
% This method computes the gradient and Hessian of a barrier function term
% corresponding to a single weight for the polynomial optimization problem.
% --------------------------------------------------------------------------
% USAGE of "gH_SOSWt"
% [in, g, H] = gH_SOSWt(x, P)
% --------------------------------------------------------------------------
% INPUT
% x:    primal iterate
% P:    evaluations of "weighted" basis polynomials at the
%       interpolation points
%
% OUTPUT
% in:   0 if P'*diag(x)*P is not positive definite. 1 if P'*diag(x)*P is
%       positive definite.
% g:    gradient of the barrier function term corresponding to a single
%       weight at x
% H:    Hessian of the barrier function term corresponding to a single
%       weight at x
% --------------------------------------------------------------------------
% EXTERNAL FUNCTIONS CALLED IN THIS FUNCTION
% None.
% -------------------------------------------------------------------------

    Y = P'*diag(x)*P;

    if ~issymmetric(Y)
        Y = (Y+Y')/2;
    end

    [L,err] = chol(Y,'lower');
    if err > 0
        in = 0;
        g = NaN;
        H = NaN;
    else
        in = 1;
        V = L\P';
        VtV = V'*V;

        g = -diag(VtV);
        H = VtV.^2;
    end
    
return
