%% -------------------------gradientDescentSolver.m--------------------------------
%
% General routine used by phase retrieval algorithms that function by using
% line search methods. This function is internal and should not be called
% by any code outside of this software package.

% The line search approach first finds a descent direction along which the
% objective function f will be reduced and then computes a step size that
% determines how far x  should move along that direction. The descent
% direction can be computed by various methods, such as gradient descent,
% Newton's method and Quasi-Newton method. The step size can be determined
% either exactly or inexactly.

% This line search algorithm implements the steepest descent, non linear
% conjugate gradient, and the LBFGS method. Set the option accordingly as
% described below.
%
%% Aditional Parameters
% The following are additional parameters that are to be passed as fields
% of the struct 'opts':
%
%     maxIters (required) - The maximum number of iterations that are
%     allowed to
%         occur.
%
%     maxTime (required) - The maximum amount of time in seconds the
%     algorithm
%         is allowed to spend solving before terminating.
%
%     tol (required) - Positive real number representing how precise the
%     final
%         estimate should be. Lower values indicate to the solver that a
%         more precise estimate should be obtained.
%
%     verbose (required) - Integer representing whether / how verbose
%         information should be displayed via output. If verbose == 0, no
%         output is displayed. If verbose == 1, output is displayed only
%         when the algorithm terminates. If verbose == 2, output is
%         displayed after every iteration.
%
%     recordTimes (required) - Whether the algorithm should store the total
%         processing time upon each iteration in a list to be obtained via
%         output.
%
%     recordResiduals (required) - Whether the algorithm should store the
%         relative residual values upon each iteration in a list to be
%         obtained via output.
%
%     recordMeasurementErrors (required) - Whether the algorithm should
%     store
%         the relative measurement errors upon each iteration in a list to
%         be obtained via output.
%
%     recordReconErrors (required) - Whether the algorithm should store the
%         relative reconstruction errors upon each iteration in a list to
%         be obtained via output. This parameter can only be set 'true'
%         when the additional parameter 'xt' is non-empty.
%
%     xt (required) - The true signal to be estimated, or an empty vector
%     if the
%         true estimate was not provided.
%
%     searchMethod (optional) - A string representing the method used to
%         determine search direction upon each iteration. Must be one of
%         {'steepestDescent', 'NCG', 'LBFGS'}. If equal to
%         'steepestDescent', then the steepest descent search method is
%         used. If equal to 'NCG', a nonlinear conjugate gradient method is
%         used. If equal to 'LBFGS', a Limited-Memory BFGS method is used.
%         Default value is 'steepestDescent'.
%
%     updateObjectivePeriod (optional) - The maximum number of iterations
%     that
%         are allowed to occur between updates to the objective function.
%         Default value is infinite (no limit is applied).
%
%     tolerancePenaltyLimit (optional) - The maximum tolerable penalty
%     caused by
%         surpassing the tolerance threshold before terminating. Default
%         value is 3.
%
%     betaChoice (optional) - A string representing the choice of the value
%         'beta' when a nonlinear conjugate gradient method is used. Must
%         be one of {'HS', 'FR', 'PR', 'DY'}. If equal to 'HS', the
%         Hestenes-Stiefel method is used. If equal to 'FR', the
%         Fletcher-Reeves method is used. If equal to 'PR', the
%         Polak-Ribi�re method is used. If equal to 'DY', the Dai-Yuan
%         method is used. This field is only used when searchMethod is set
%         to 'NCG'. Default value is 'HS'.
%
%     ncgResetPeriod (optional) - The maximum number of iterations that are
%         allowed to occur between resettings of a nonlinear conjugate
%         gradient search direction. This field is only used when
%         searchMethod is set to 'NCG'. Default value is 100.
%
%     storedVectors (optional) - The maximum number of previous iterations
%     of
%         which to retain LBFGS-specific iteration data. This field is only
%         used when searchMethod is set to 'LBFGS'. Default value is 5.

function [sol, outs] = gradientDescentSolver(A, At, x0, b0, updateObjective, opts)
setDefaultOpts();

% Length of input signal
n = length(x0);
if ~isempty(opts.xt)
    residualTolerance = 1.0e-13;
else
    residualTolerance = opts.tol;
end

% Iteration number of last objective update
lastObjectiveUpdateIter = 0;
% Total penalty caused by surpassing tolerance threshold
tolerancePenalty = 0;
% Whether to update objective function upon next iteration
updateObjectiveNow = true;
% Maximum norm of differences between consecutive estimates
maxDiff = -inf;

currentSolveTime = 0;
currentMeasurementError = [];
currentResidual = [];
currentReconError = [];

if opts.recordTimes
    solveTimes = zeros(opts.maxIters, 1);
end
if opts.recordResiduals
    residuals = zeros(opts.maxIters, 1);
end
if opts.recordMeasurementErrors
    measurementErrors = zeros(opts.maxIters, 1);
end
if opts.recordReconErrors
    reconErrors = zeros(opts.maxIters, 1);
end

x1 = x0;
d1 = A(x1);

startTime = tic;
for iter = 1 : opts.maxIters
    % Signal to update objective function after fixed number of iterations
    % have passed
    if iter - lastObjectiveUpdateIter == opts.updateObjectivePeriod
        updateObjectiveNow = true;
    end
    % Update objective if flag is set
    if updateObjectiveNow
        updateObjectiveNow = false;
        lastObjectiveUpdateIter = iter;
        [f, gradf] = updateObjective(x1, d1);
        f1 = f(d1);
        gradf1 = At(gradf(d1));
        
        if strcmpi(opts.searchMethod, 'lbfgs')
            % Perform LBFGS initialization
            yVals = zeros(n, opts.storedVectors);
            sVals = zeros(n, opts.storedVectors);
            rhoVals = zeros(1, opts.storedVectors);
        elseif strcmpi(opts.searchMethod, 'ncg')
            % Perform NCG initialization
            lastNcgResetIter = iter;
            unscaledSearchDir = zeros(n, 1);
        end
        searchDir1 = determineSearchDirection();
        % Reinitialize stepsize to supplement new objective function
        tau1 = determineInitialStepsize();
    else
        gradf1 = At(gradf(d1));
        Dg = gradf1 - gradf0;
        
        if strcmpi(opts.searchMethod, 'lbfgs')
            % Update LBFGS stored vectors
            sVals = [Dx, sVals(:, 1:opts.storedVectors-1)];
            yVals = [Dg, yVals(:, 1:opts.storedVectors-1)];
            rhoVals = [1 / real(Dg' * Dx), rhoVals(:, 1:opts.storedVectors-1)];
        end
        
        searchDir1 = determineSearchDirection();
        updateStepsize();
    end
    
    x0 = x1;
    f0 = f1;
    gradf0 = gradf1;
    tau0 = tau1;
    searchDir0 = searchDir1;
    
    x1 = x0 + tau0 * searchDir0;
    Dx = x1 - x0;
    d1 = A(x1);
    f1 = f(d1);
    
    % We now determine an appropriate stepsize for our algorithm using
    % Armijo-Goldstein condition
    
    backtrackCount = 0;
    % Cap maximum number of backtracks
    while backtrackCount <= 20
        tmp = f0 + 0.1 * tau0 * real(searchDir0' * gradf0);
        
        % Break if f1 < tmp or f1 is sufficiently close to tmp (determined
        % by error)
        % Avoids division by zero
        if f1 <= tmp 
            break;
        end
        
        backtrackCount = backtrackCount + 1;
        % Stepsize reduced by factor of 5
        tau0 = tau0 * 0.2;
        x1 = x0 + tau0 * searchDir0;
        Dx = x1 - x0;
        d1 = A(x1);
        f1 = f(d1);
    end
    
    % Handle processing of current iteration estimate
    stopNow = processIteration();
    if stopNow
        break;
    end
end

% Create output
sol = x1;
outs = struct;
outs.iterationCount = iter;
if opts.recordTimes
    outs.solveTimes = solveTimes;
end
if opts.recordResiduals
    outs.residuals = residuals;
end
if opts.recordMeasurementErrors
    outs.measurementErrors = measurementErrors;
end
if opts.recordReconErrors
    outs.reconErrors = reconErrors;
end

if opts.verbose == 1
    displayVerboseOutput();
end

% Assigns default options for any options that were not provided by the
% client
    function setDefaultOpts()
        if ~isfield(opts, 'updateObjectivePeriod')
            % Objective function is never updated by default
            opts.updateObjectivePeriod = inf;
        end
        if ~isfield(opts, 'tolerancePenaltyLimit')
            opts.tolerancePenaltyLimit = 3;
        end
        if ~isfield(opts, 'searchMethod')
            opts.searchMethod = 'steepestDescent';
        end
        if strcmpi(opts.searchMethod, 'lbfgs')
            if ~isfield(opts, 'storedVectors')
                opts.storedVectors = 5;
            end
        end
        if strcmpi(opts.searchMethod, 'ncg')
            if ~isfield(opts, 'betaChoice')
                opts.betaChoice = 'HS';
            end
            if ~isfield(opts, 'ncgResetPeriod')
                opts.ncgResetPeriod = 100;
            end
        end
    end

% Determine reasonable initial stepsize of current objective function
% (adapted from FASTA.m)
    function tau = determineInitialStepsize()
        x_1 = randn(size(x0));
        x_2 = randn(size(x0));
        gradf_1 = At(gradf(A(x_1)));
        gradf_2 = At(gradf(A(x_2)));
        L = norm(gradf_1-gradf_2) / norm(x_2-x_1);
        L = max(L, 1.0e-30);
        tau = 25.0 / L;
    end

% Determine search direction for next iteration based on specified search
% method
    function searchDir = determineSearchDirection()
        switch lower(opts.searchMethod)
            case 'steepestdescent'
                searchDir = -gradf1;
            case 'ncg'
                searchDir = -gradf1;
                
                % Reset NCG progress after specified number of iterations have
                % passed
                if iter - lastNcgResetIter == opts.ncgResetPeriod
                    unscaledSearchDir = zeros(n, 1);
                    lastNcgResetIter = iter;
                end
                
                % Proceed only if reset has not just occurred
                if iter ~= lastNcgResetIter
                    switch lower(opts.betaChoice)
                        case 'hs'
                            % Hestenes-Stiefel
                            beta = -real(gradf1' * Dg) / real(unscaledSearchDir' * Dg);
                        case 'fr'
                            % Fletcher-Reeves
                            beta = norm(gradf1)^2 / norm(gradf0)^2;
                        case 'pr'
                            % Polak-Ribi�re
                            beta = real(gradf1' * Dg) / norm(gradf0)^2;
                        case 'dy'
                            % Dai-Yuan
                            beta = norm(gradf1)^2 / real(unscaledSearchDir' * Dg);
                    end
                    searchDir = searchDir + beta * unscaledSearchDir;
                end
                unscaledSearchDir = searchDir;
            case 'lbfgs'
                searchDir = -gradf1;
                iters = min(iter - lastObjectiveUpdateIter, opts.storedVectors);
                
                if iters > 0
                    alphas = zeros(iters, 1);
                    
                    % First loop
                    for j = 1 : iters
                        alphas(j) = rhoVals(j) * real(sVals(:,j)' * searchDir);
                        searchDir = searchDir - alphas(j) * yVals(:,j);
                    end
                    
                    % Scaling of search direction
                    gamma = real(Dg' * Dx) / (Dg' * Dg);
                    searchDir = gamma * searchDir;
                    
                    % Second loop
                    for j = iters : -1 : 1
                        beta = rhoVals(j) * real(yVals(:,j)' * searchDir);
                        searchDir = searchDir + (alphas(j) - beta) * sVals(:,j);
                    end
                    
                    searchDir = 1/gamma * searchDir;
                    searchDir = norm(gradf1) / norm(searchDir) * searchDir;
                end
        end
        
        % Change search direction to steepest descent direction if current
        % direction is invalid
        if any(isnan(searchDir)) || any(isinf(searchDir))
            searchDir = -gradf1;
        end
        
        % Scale current search direction match magnitude of gradient
        searchDir = norm(gradf1) / norm(searchDir) * searchDir;
    end

% Update stepsize when objective update has not just occurred (adopted from
% FASTA.m)
    function updateStepsize()
        Ds = searchDir0 - searchDir1;
        dotprod = real(dot(Dx, Ds));
        tauS = norm(Dx)^2 / dotprod;  %  First BB stepsize rule
        tauM = dotprod / norm(Ds)^2; %  Alternate BB stepsize rule
        tauM = max(tauM, 0);
        if 2*tauM > tauS   %  Use "Adaptive"  combination of tau_s and tau_m
            tau1 = tauM;
        else
            tau1 = tauS - tauM / 2;  %  Experiment with this param
        end
        if tau1 <= 0 || isinf(tau1) || isnan(tau1)      %  Make sure step is non-negative
            tau1 = tau0 * 1.5;  % let tau grow, backtracking will kick in if stepsize is too big
        end
    end

    function stopNow = processIteration()
        currentSolveTime = toc(startTime);
        maxDiff = max(norm(Dx), maxDiff);
        currentResidual = norm(Dx) / maxDiff;
        
        % Obtain recon error only if true solution has been provided
        if ~isempty(opts.xt)
            reconEstimate = (x1'*opts.xt)/(x1'*x1) * x1;
            currentReconError = norm(opts.xt - reconEstimate) / norm(opts.xt);
        end
        
        if opts.recordTimes
            solveTimes(iter) = currentSolveTime;
        end
        if opts.recordResiduals
            residuals(iter) = currentResidual;
        end
        if opts.recordMeasurementErrors
            currentMeasurementError = norm(abs(d1) - b0) / norm(b0);
            measurementErrors(iter) = currentMeasurementError;
        end
        if opts.recordReconErrors
            assert(~isempty(opts.xt),['You must specify the ground truth solution '...
                ,'if the "recordReconErrors" flag is set to true.  Turn '...
                ,'this flag off, or specify the ground truth solution.']);
            reconErrors(iter) = currentReconError;
        end
        
        if opts.verbose == 2
            displayVerboseOutput();
        end
        
        % Terminate if solver surpasses maximum allocated timespan
        if currentSolveTime > opts.maxTime
            stopNow = true;
            return;
        end
        
        % If user has supplied actual solution, use recon error to determine
        % termination
        if ~isempty(currentReconError)
            if currentReconError <= opts.tol
                stopNow = true;
                return;
            end
        end

        if currentResidual <= residualTolerance
            % Give algorithm chance to recover if stuck at local minimum by
            % forcing update of objective function
            updateObjectiveNow = true;
            % If algorithm continues to get stuck, terminate
            tolerancePenalty = tolerancePenalty + 1;
            if tolerancePenalty >= opts.tolerancePenaltyLimit
                stopNow = true;
                return;
            end
        end
        
        stopNow = false;
    end

% Display output to user based on provided options
    function displayVerboseOutput()
        fprintf('Iter = %d', iter);
        fprintf(' | IterationTime = %.3f', currentSolveTime);
        fprintf(' | Resid = %.3e', currentResidual);
        fprintf(' | Stepsize = %.3e', tau0);
        if ~isempty(currentMeasurementError)
            fprintf(' | M error = %.3e', currentMeasurementError);
        end
        if ~isempty(currentReconError)
            fprintf(' | R error = %.3e', currentReconError);
        end
        fprintf('\n');
    end
end