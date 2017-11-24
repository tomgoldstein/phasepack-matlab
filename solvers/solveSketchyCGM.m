%                           solveSketchyCGM.m
%
%  Solver for the skechyCGM algorithm.
%  
%% I/O
%  Inputs:
%     A:    m x n matrix or a function handle to a method that
%           returns A*x.
%     At:   The adjoint (transpose) of 'A'. If 'A' is a function handle,
%           'At' must be provided.
%     b0:   m x 1 real,non-negative vector consists of all the
%           measurements. 
%     x0:   n x 1 vector. It is the initial guess of the unknown signal x.
%     opts: A struct consists of the options for the algorithm. 

% For details, see header in solvePhaseRetrieval.m or the User Guide.
%
%     Note: When a function handle is used, the value of 'At' (a function
%     handle for the adjoint of 'A') must be supplied.
% 
%  Outputs:
%     sol:  n x 1 vector. It is the estimated signal. 
%     outs: A struct consists of the convergence info. 

% For details, see header in solvePhaseRetrieval.m or the User Guide.
%  
%  
%  See the script 'testSketchyCGM.m' for an example of proper usage of this
%  function.
%
%% Notations
%  B.3.1 X = xx' is the lifted version of the unknown signal x. The linear
%  operator curly_A: R(nxn) -> R(d). And its adjoint, curly_As: R(d) ->
%  R(nxn)
%                 curly_A(X) = diag(A*X*A') curly_As(z) = A' * diag(z)' A
%
%
%% Algorithm Description
%  SketchyCGM modifies a standard convex optimization scheme, the
%  conditional gradient method, to store only a small randomized sketch of
%  the matrix variable. After the optimization terminates, the algorithm
%  extracts a low-rank approximation of the solution from the sketch.

%  This algorithm solves the problem:
%                 minimize f(curly_A(X)) s.t. ||X||_N <= \alpha where
%                 ||X||_N is the nuclear norm of X i.e. trace(X) and \alpha
%                 is a constant value (see section 1.2 of SketchCGM paper
%                 for it).
%  
%  Specifically, the method has the following five steps(section 5.1 of the
%  paper): (1) Initialize iterate z, two random matrices Psi and Omega, two
%  sketches Y and W. (2) At each iteration, compute an update direction via
%  Lanczos or via randomized SVD (3) Update the iterate and the two
%  sketches Y and W (4) The iteration continues until it triggers the
%  stopping condition(i.e. the normalized gradient is smaller than a given
%  tolerance). (5) Form a rank-r approximate solution X of the model
%  problem and reconstruct final solution x from it.
%  
%  Note: The suboptimality eps used in the paper for stopping condition
%  doesn't work well in our tests so we use normalized gradient instead.
%  
%  For a detailed explanation, see the paper referenced below.
%

%% References
%  Paper Title:   Sketchy Decisions: Convex Low-Rank Matrix Optimization
%  with Optimal Storage Place:         Chapter 5.4, Algorithm 1 Authors:
%  Alp Yurtsever, Madeleine Udell, Joel A. Tropp, Volkan Cevher arXiv
%  Address: https://arxiv.org/abs/1702.06838
%  
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein Copyright (c) University of Maryland,
% 2017


%% -----------------------------START----------------------------------- 
 

function [sol,outs] = solveSketchyCGM(A,At,b0,x0,opts)

    validateOptions(opts); % Check the validify of algorithm specific options

    % Initialize test setup
    n = length(x0);  % length of unknown signal
    m = length(b0);  % number of measurements
    b = b0.^2;       % to be consistent with the paper
    r = opts.rank;   % the rank parameter that influences the dimensions of the 
                     % "sketch" Y and W and random matrices Psi and Omega. 
                     % So the approximate solution will have rank r.
    z = zeros(m, 1); % Initialize the iterate
    alpha = sum(b)/m;% See section 1.2 of SketchCGM paper
    maxNormGradfz = 1e-30; % Initialize the max norm of gradient to be a small non-zero value 

    % Initialize values potentially computed at each round.
    currentTime = [];
    currentResid = [];
    currentReconError = [];
    currentMeasurementError = [];
    sol = x0;
 
    % Initialize vectors for recording convergence information
    [solveTimes,measurementErrors,reconErrors,residuals] = initializeContainers(opts);

    % Initialize Sketchy parameters
    [Omega, Psi, Y, W] = sketchyInit(m, n, r);

    % Build linear operators
    curly_A_of_uut = @(u) sum(abs(A(u)).^2,2); % Compute curly_A_of_uut(u)=curly_A(X)=diag(A*X*A')
                                               % where X = uu'
    curly_As = @(z) @(x) At(z.*A(x));          % Given vector 'z', return function handle
                                               % for the operator/matrix A'*diag(z)*A
    % gradient of the objective function
    % .5 * ||z - b||^2
    grad_f = @(z) z - b;                       
    
    % record starting point
    Y0=Y; W0=W;
        
    % stepsize
    eta = opts.eta;
    failCounter = 0;
    startTime = tic; % Start timer

    % Start SketchyCGM iterations
    for iter = 1 : opts.maxIters
        grad_fz = grad_f(z);
        normGradfz = norm(grad_fz);
        if normGradfz > maxNormGradfz % Update max norm of gradient if necessary
            maxNormGradfz = normGradfz;
        end

        % Return the minimal eigenvalue and eigenvector
        [lambda,u] = MinEig(curly_As(grad_fz),n);
        
        % Update Rule (from Section 5.5 of SketchCGM paper)
        if lambda>0
            h=0;  % In the SketchyCGM paper, the authors say to set this to zero
            u=0;  % You ALSO need to set this to zero, but the authors neglect to mention this
        else
            h = alpha*curly_A_of_uut(u);
        end
            

        % Else continue update
        z = (1 - eta)*z + eta*h;
        grad_fz = grad_f(z);
        
        % Update the sketches
        [Y, W] = sketchyUpdate(-alpha*u, -u, eta, Y, W, Omega, Psi);

        % Calculate residual
        currentResid = norm(grad_fz)/maxNormGradfz;
        residuals(iter) = currentResid;
        % Cut stepsize if residual goes up
        
        
        % Select stepsize
        % During optimization, cut the stepsize if the errors don't go
        % down.
        if iter>1 && residuals(iter)>residuals(iter-1)
            failCounter = failCounter+1;
        end
        if failCounter>3
            eta = eta/2;
            failCounter = 0;
        end
        
        
        %---------------------------------------------------------------------------
        % Record convergence information and check stopping condition
        % If xt is provided, reconstruction error will be computed and used for stopping
        % condition. Otherwise, residual will be computed and used for stopping
        % condition.
        if ~isempty(opts.xt)
            if iter == 1
                sol = x0;
            else
                sol = recoverSignal(Y, W, Psi, r);
            end
                        
            x = sol;
            xt = opts.xt;
            %  Compute optimal rotation
            alpha2 = (x(:)'*xt(:))/(x(:)'*x(:));
            x = alpha2*x;
            currentReconError = norm(x-xt)/norm(xt);
            if opts.recordReconErrors
                reconErrors(iter) = currentReconError;
            end
        end

        currentTime = toc(startTime);                % Measure elapsed time so far
        
        if opts.recordTimes
            solveTimes(iter) = currentTime;
        end

        if opts.recordMeasurementErrors
            if iter == 1
                sol = x0;
            else
                sol = recoverSignal(Y, W, Psi, r);
            end
            currentMeasurementError = norm(abs(A(sol)) - b0) / norm(b0);
            measurementErrors(iter) = currentMeasurementError;
        end
       
        % Display verbose output if specified
        if opts.verbose == 2
            displayVerboseOutput(iter, currentTime, currentResid, currentReconError, currentMeasurementError);
        end

        %  Test stopping criteria. 
        if stopNow(opts, currentTime, currentResid, currentReconError)
            break;
        end
        %---------------------------------------------------------------------------

    end
    

    % Recover the signal from two sketches Y, W and matrix Psi
    sol = recoverSignal(Y, W, Psi, r);
    % Create output according to the options chosen by user
    outs = generateOutputs(opts, iter, solveTimes, measurementErrors, reconErrors, residuals);

    % Display verbose output if specified
    if opts.verbose == 1
        displayVerboseOutput(iter, currentTime, currentResid, currentReconError, currentMeasurementError);
    end
end


%% Helper Functions
% Initialize Sketchy parameters: random matrices Omega and Psi, sketches Y and W
function [Omega, Psi, Y, W] = sketchyInit(m, n, r)
    k = 2*r + 1;
    l = 4*r + 3;
    Omega = randn(n,k);
    Psi = randn(l,n);
    Y = zeros(n,k);
    W = zeros(l,n);
end

% Update the sketch
function [Y, W] = sketchyUpdate(u, v, eta, Y, W, Omega, Psi)
    Y = (1 - eta)*Y + eta*u*(v'*Omega);
    W = (1 - eta)*W + eta*(Psi*u)*v';
end

% Reconstruct the lifted solution U*S*V = xx' from sketch
function [U, S, V] = sketchyReconstruct(Y, W, Psi, r)
    Q = orth(Y);
    B = (Psi*Q)\W;
    [U, S, V] = svds(B, r);
    U = Q*U;
end


%  Return the minimal (hopefully negative) eigenvalue and eigenvector
function [lambda,u] = MinEig(matrix,n)
    eigs_opts.isreal = false;
    [u, lambda] = eigs(matrix,n,1,'LM',eigs_opts); %  get largest magnitude eigenvector
    if lambda>0    % check that largest magnitude vector was also the minimum vector
        [u, lambda1] = eigs(@(x) matrix(x) - lambda*x,n,1,'LM',eigs_opts);
        lambda = lambda1+lambda;
    end
end 

% Recover the signal from two sketches Y, W and matrix Psi
function sol = recoverSignal(Y, W, Psi, r)
    % Reconstruct the lifted solution X_hat = xx*
    [U, S, V] = sketchyReconstruct(Y, W, Psi, r);
    % Reconstruct final solution 
    sol = U*sqrt(S);
end

% Check the validify of algorithm specific options
function validateOptions(opts)
    checkIfNumber('rank for sketchuyCGM', opts.rank);
end
