%                           solvePhaseLift.m
%
%  Implementation of the PhaseLift algorithm using full-scale semidefinite
%  programming.
%
%% I/O
%  Inputs:
%     A:    m x n matrix or a function handle to a method that
%           returns A*x.     
%     At:   The adjoint (transpose) of 'A'. If 'A' is a function handle,
%           'At' must be provided.
%     b0:   m x 1 real,non-negative vector consists of all the measurements.
%     x0:   n x 1 vector. It is the initial guess of the unknown signal x.
%     opts: A struct consists of the options for the algorithm. For details,
%           see header in solvePhaseRetrieval.m or the User Guide.
%
%     Note: When a function handle is used, the
%     value of 'At' (a function handle for the adjoint of 'A') must be 
%     supplied.
% 
%  Outputs:
%     sol:  n x 1 vector. It is the estimated signal.
%     outs: A struct consists of the convergence info. For details,
%           see header in solvePhaseRetrieval.m or the User Guide.
%  
%  
%  See the script 'runPhaseLift.m' for an example of proper usage of this 
%  function.
%
%% Notations
%  X = x*x' is the lifted version of the unknown signal x.
%  b = b0.^2 is the element-wise square of the measurements b0.
%  AL is a function handle that takes X as input and outputs b.
%  
%% Algorithm Description
%  PhaseLift algorithm reformualates the PR problem as a convex problem by
%  lifting the dimension of the unknown signal x. The problem becomes:
%  minimize Tr(X) subject to AL(X)=b and X is a positive semidefinite
%  matrix.
%  
%  More specifically,
%  Solve the problem
%           min  mu||X||_nuc +.5||A(X)-b||^2 
%                X>=0
%  where X is a square symmetric matrix,||X||_nuc is the nuclear (trace) 
%  norm, A is a linear operator, and X>=0 denotes that X
%  must lie in the positive semidefinite cone.
%
%  The unknown signal x can be recovered from its lifted version X by
%  factorization.
%
%  The solver has three steps
%  (1) Lifting x to X=xx' and create AL and its transpose AtL using A and At.
%  (2) Use FASTA to solve this convex optimization problem.
%  (3) Take the principle eigenvector of X and re-scale it to capture 
%      the correct ammount of energy.
%  
%  For a detailed explanation, see the PhaseLift paper referenced below.
%  For more details about FASTA, see the FASTA user guide, or the paper "A
%  field guide to forward-backward splitting with a FASTA implementation."
%
%% References
%  Paper Title:   PhaseLift: Exact and Stable Signal Recovery from
%                 Magnitude Measurements via Convex Programming
%  Place:         Chapter 2.3
%  Authors:       Emmanuel J. Candes, Thomas Strohmer, Vladislav Voroninski
%  arXiv Address: https://arxiv.org/abs/1109.4499
%  
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017




%% -----------------------------START----------------------------------- 
 

function [sol,outs] = solvePhaseLift(A, At, b0, x0, opts)

    % Initialization
    muFinal = opts.regularizationPara;   % Regularization Parameter
    m = length(b0);                 % Number of unknowns
    n = numel(x0);                  % The dimension of the lifted matrix
    b = b0.^2;                      % To be consistent with the notation in the paper
    X0 = x0*x0';                    % Lifted Initial Guess

    %  We need direct access to the entires of A, so convert a function
    %  handle to a dense matrix
    if ~isnumeric(A)
        A = A(eye(size(x0,1)));
    end
   
    % Initialize values potentially computed at each round.
    iter = 0;
    currentTime = [];
    currentResid = [];
    currentReconError = [];
    currentMeasurementError = [];

    % Initialize vectors for recording convergence information
    [solveTimes,measurementErrors,reconErrors,residuals] = initializeContainers(opts);

    % Create wrappers for functions passed to FASTA
    % These function compute the linear measurement operator on the lifted
    % matrix, and it's adjoint
    AL = @(X) sum((A*X).*conj(A),2);
    AtL = @(y)  A'*(A.*(y*ones(1,n)));
    
    % Choose a regulatization parameter (the coefficient of the nuclear
    % norm term).  Note that norm(eig(AtL(b0)),1) is the regularizer that
    % produces a zero solution.
    mu = norm(eig(AtL(b0)),1)*opts.tol;
 
    
     % The current loop function handle version is only about 10% slower 
    % than the vectorized matrix version.
    % It can probably be further optimized by
    % 1.loop roll-out.
    % 2.AL0 can take advantage of the factorization of X.
    
    %%  Define ingredients for FASTA
    % m x 1 -> 1 x 1
    % f(y) = .5 ||y - b||^2
    f = @(y) .5*norm(y(:)-b(:))^2;

    % m x 1 -> m x 1
    gradf = @(y) y(:)-b(:);

    % n x n -> 1 x 1
    % g(z) = mu||X||_nuc, plus characteristic function of the SDP cone
    g  = @(X) mu*norm(eig(X),1);
    
    % n x n -> n x n
    % proxg(z,t) = argmin t*mu*nuc(x)+.5||x-z||^2, with x in SDP cone
    proxg = @(X,t)  projectSemiDefCone(X, mu*t);

    % Options to pass to fasta
    fastaOpts.tol = opts.tol;
    fastaOpts.maxIters = opts.maxIters;
    fastaOpts.accelerate = false;
    fastaOpts.stopNow = @(x, iter, resid, normResid, maxResid, opts) ...
        processIteration(x, normResid);  % Use customized stopNow in order to get
                                         % solveTime, residual and error at each iteration.


    %% Call solver
    startTime = tic;
    [Xest, outs] = fasta(AL,AtL,f,gradf,g,proxg,X0,fastaOpts);

    sol = recoverSignal(Xest,n);

    % Create output according to the options chosen by user
    outs = generateOutputs(opts, iter, solveTimes, measurementErrors, reconErrors, residuals);

    % Display verbose output if specified
    if opts.verbose == 1
        displayVerboseOutput(iter, currentTime, currentResid, currentReconError, currentMeasurementError);
    end


    % Runs code upon each FASTA iteration. Returns whether FASTA should
    % terminate.
    function stop = processIteration(x, normResid)
        iter = iter + 1;
        % Record convergence information and check stopping condition
        % If xt is provided, reconstruction error will be computed and used for stopping
        % condition. Otherwise, residual will be computed and used for stopping
        % condition.
        if ~isempty(opts.xt)
            
            xt = opts.xt;
            xt = xt*xt'; 
            %  Compute optimal rotation
            alpha = (x(:)'*xt(:))/(x(:)'*x(:));
            x = alpha*x;
            currentReconError = norm(x-xt)/norm(xt);
            if opts.recordReconErrors
                reconErrors(iter) = currentReconError;
            end
        end

        if isempty(opts.xt) | opts.recordResiduals
            currentResid = normResid;
        end

        if opts.recordResiduals
            residuals(iter) = currentResid;
        end

        currentTime = toc(startTime);                % Record elapsed time so far
        if opts.recordTimes
            solveTimes(iter) = currentTime;
        end
        if opts.recordMeasurementErrors
            currentMeasurementError = norm(abs(A(sol)) - b0) / norm(b0);
            measurementErrors(iter) = currentMeasurementError;
        end
       
        % Display verbose output if specified
        if opts.verbose == 2
            displayVerboseOutput(iter, currentTime, currentResid, currentReconError, currentMeasurementError);
        end

        % Test stopping criteria. 
        stop = stopNow(opts, currentTime, currentResid, currentReconError);
    end

    %% Recover solution using method recommended by PhaseLift authors
    %  The solution matrix may not be rank 1.  In this case, we use the
    %  principle eigenvector and re-scale it to capture the correct ammount 
    %  of energy.
    function [sol]=recoverSignal(Xest, n)
        %  Extract X
        [V,D] = eig(reshape(Xest,[n,n])); %  get principle eigenvector
        [val,ind] = max(diag(D));
        recovered = V(:,ind)*sqrt(D(ind,ind));
        % figure out how muhc energy we're missing because the solution might
        % have other eigenvectors
        lifted = recovered*recovered';
        scale = norm(b)/norm(AL(lifted));
        % scale the solution to capture the lost energy
        sol = recovered*scale;
    end
end
% n x n -> n x n
% proxg(z,t) = argmin t*mu*nuc(x)+.5||x-z||^2, with x in SDP cone
function [X] = projectSemiDefCone(X, delta)
    [V,D] = eig(X);
    D = max(real(D)-delta,0);
    X = V*D*V';
end

% Check the validify of algorithm specific options
function validateOptions(opts)
    checkIfNumber('regularizationPara', opts.regularizationPara);
end

