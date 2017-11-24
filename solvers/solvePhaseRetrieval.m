%                               solvePhaseRetrieval.m
%   This method solves the problem:

%                        Find x given b0 = |Ax+epsilon|
%   
%  Where A is a m by n complex matrix, x is a n by 1 complex vector, b0 is
%   a m by 1 real,non-negative vector and epsilon is a m by 1 vector. The
%   user supplies function handles A, At and measurement b0. Note: The
%   unknown signal to be recovered must be 1D for our interface. Inputs:
%    A     : A m x n matrix (or optionally a function handle to a method)
%            that returns A*x
%    At    : The adjoint (transpose) of 'A.' It can be a n x m matrix or a
%            function handle. 
%    b0    : A m x 1 real,non-negative vector consists of  all the
%            measurements. 
%    n     : The size of the unknown signal. It must be provided if A is a
%            function handle.
%    opts  : An optional struct with options.  The commonly used fields
%            of 'opts' are:
%               initMethod              : (string,
%               default='truncatedSpectral') The method used
%                                         to generate the initial guess x0.
%                                         User can use a customized initial
%                                         guess x0 by providing value to
%                                         the field customx0.
%               algorithm               : (string, default='altmin') The
%                                         algorithm used
%                                         to solve the phase retrieval
%                                         algorithm. User can use a
%                                         customized algorithm by providing
%                                         a function [A,At,b0,x0,opts]
%                                         ->[x,outs,opts] to the field
%                                         customAlgorithm.
%               maxIters                : (integer, default=1000) The
%                                         maximum number of
%                                         iterations allowed before
%                                         termination.
%               maxTime                 : (positive real number,
%                                         default=120, unit=second)
%                                         The maximum seconds allowed
%                                         before termination.
%               tol                     : (double, default=1e-6) The
%                                         stopping tolerance.
%                                         It will be compared with
%                                         reconerror if xt is provided.
%                                         Otherwise, it will be compared
%                                         with residual. A smaller value of
%                                         'tol' usually results in more
%                                         iterations.
%               verbose                 : ([0,1,2], default=0)  If ==1,
%                                         print out
%                                         convergence information in the
%                                         end. If ==2, print out
%                                         convergence information at every
%                                         iteration.
%               recordMeasurementErrors : (boolean, default=false) Whether
%                                         to compute and record
%                                         error(i.e.
%                                         norm(|Ax|-b0)/norm(b0)) at each
%                                         iteration.
%               recordResiduals         : (boolean, default=true) If it's
%                                         true, residual will be
%                                         computed and recorded at each
%                                         iteration. If it's false,
%                                         residual won't be recorded.
%                                         Residual also won't be computed
%                                         if xt is provided. Note: The
%                                         error metric of residual varies
%                                         across solvers.
%               recordReconErrors       : (boolean, default=false) Whether
%                                         to record
%                                         reconstruction error. If it's
%                                         true, opts.xt must be provided.
%                                         If xt is provided reconstruction
%                                         error will be computed regardless
%                                         of this flag and used for
%                                         stopping condition.
%               recordTimes             : (boolean, default=true) Whether
%                                         to record
%                                         time at each iteration. Time will
%                                         be measured regardless of this
%                                         flag.
%               xt                      : A n x 1 vector. It is the true
%                                         signal. Its purpose is
%                                         to compute reconerror.
%
%            There are other more algorithms specific options not listed
%            here. To use these options, set the corresponding field in
%            'opts'. For example:
%                      >> opts.tol=1e-8; >> opts.maxIters = 100;
%
%
%   Outputs:
%    sol               : The approximate solution outs : A struct with
%                        convergence information
%    iterationCount    : An integer that is  the number of iterations the
%                        algorithm runs.
%    solveTimes        : A vector consists  of elapsed (exist when
%                        recordTimes==true) time at each iteration. 
%    measurementErrors : A vector consists of the errors (exist when
%                        recordMeasurementErrors==true)   i.e.
%                        norm(abs(A*x-b0))/norm(b0) at each iteration.
%    reconErrors       : A vector consists of the reconstruction (exist
%                        when recordReconErrors==true) errors 
%                        i.e. norm(xt-x)/norm(xt) at each iteration.
%    residuals         : A vector consists of values that (exist when
%                        recordResiduals==true)  will be compared with 
%                        opts.tol for stopping condition  checking.
%                        Definition varies across solvers.
%    opts              : A struct that contains fields used by the solver. 
%                        Its possible fields are the same as the input 
%                        parameter opts.
%   
% For more details and more options in opts, see the Phasepack user
%   guide.

% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein Copyright (c) University of Maryland,
% 2017




%% -----------------------------START----------------------------------- 
 
function [sol, outs, opts] = solvePhaseRetrieval(A, At, b0, n, opts)
% Add path to helper functions
% addpath('util');
% addpath('initializers');

% If opts is not provided, create it
if ~exist('opts', 'var')
    opts = struct;
end

% If A is a matrix, infer n and At from A
if isnumeric(A) & ~isempty(A)
    n = size(A, 2);
    % Transform matrix into function form
    At = @(x) A' * x;
    A = @(x) A * x;
end

% Apply default options and validate user-defined options
opts = manageOptions(opts);
% Check that inputs are of valid datatypes and sizes
validateInput(A, At, b0, n, opts);
% Check that At is the adjoint/transpose of A
checkAdjoint(A, At, b0);

 
x0 = initX(A, At, b0, n, opts); % Initialize x0
% Truncate imaginary components of x0 if working with real values
if ~opts.isComplex
    x0 = real(x0);
else
    if opts.isNonNegativeOnly
        warning('opts.isNonNegativeOnly will not be used when the signal is complex.');
    end
end


[sol, outs] = solveX(A, At, b0, x0, opts); % Solve the problem using the specified algorithm
end

%% Helper functions

% Initialize x0 using the specified initialization method
function x0 = initX(A, At, b0, n, opts)
switch lower(opts.initMethod)
    case {'truncatedspectral','truncated'}
        x0 = initSpectral(A,At,b0,n,true,true,opts.verbose);  % isTruncated, isScaled?
    case 'spectral'
        x0 = initSpectral(A,At,b0,n,false,true,opts.verbose); % isTruncated, isScaled?
    case {'amplitudespectral','amplitude'}
        x0 = initAmplitude(A,At,b0,n,opts.verbose);
    case {'weightedspectral','weighted'}
        x0 = initWeighted(A,At,b0,n,opts.verbose);
    case {'orthogonalspectral','orthogonal'}
        x0 = initOrthogonal(A,At,b0,n,opts.verbose);
    case {'optimal','optimalspectral'}
        x0 = initOptimalSpectral(A,At,b0,n,true,opts.verbose);
    case 'angle'
        assert(isfield(opts,'xt'),'The true solution, opts.xt, must be specified to use the angle initializer.')
        assert(isfield(opts,'initAngle'),'An angle, opts.initAngle, must be specified (in radians) to use the angle initializer.')
        x0 = initAngle(opts.xt, opts.initAngle);
    case 'custom'
        x0 = opts.customx0;
    otherwise
        error('Unknown initialization method "%s"', opts.initMethod);
end
end

% Estimate x0 using the specified algorithm
function [sol, outs] = solveX(A, At, b0, x0, opts)
switch lower(opts.algorithm)
    case 'custom'
        [sol, outs] = optsCustomAlgorithm(A, At, b0, x0, opts);
    case 'amplitudeflow'
        [sol, outs] = solveAmplitudeFlow(A, At, b0, x0, opts);
    case 'coordinatedescent'
        [sol, outs] = solveCoordinateDescent(A, At, b0, x0, opts);
    case 'fienup'
        [sol, outs] = solveFienup(A, At, b0, x0, opts);
    case 'gerchbergsaxton'
        [sol, outs] = solveGerchbergSaxton(A, At, b0, x0, opts);
    case 'kaczmarz'
        [sol, outs] = solveKaczmarzSimple(A, At, b0, x0, opts);
    case 'phasemax'
        [sol, outs] = solvePhaseMax(A, At, b0, x0, opts);
    case 'phaselamp'
        [sol, outs] = solvePhaseLamp(A, At, b0, x0, opts);
    case 'phaselift'
        [sol, outs] = solvePhaseLift(A, At, b0, x0, opts);
    case 'raf'
        [sol, outs] = solveRAF(A, At, b0, x0, opts);
    case 'rwf'
        [sol, outs] = solveRWF(A, At, b0, x0, opts);
    case 'sketchycgm'
        [sol, outs] = solveSketchyCGM(A, At, b0, x0, opts);
    case 'taf'
        [sol, outs] = solveTAF(A, At, b0, x0, opts);
    case 'twf'
        [sol, outs] = solveTWF(A, At, b0, x0, opts);
    case 'wirtflow'
        [sol, outs] = solveWirtFlow(A, At, b0, x0, opts);
    otherwise
        error('Unknown algorithm "%s"', opts.algorithm);
end
end

% Check validity of input
function validateInput(A, At, b0, n, opts)
if ~isnumeric(A) & (isempty(At)|isempty(n))
    error('If A is a function handle, then At and n must be provided')
end

assert(n > 0, 'n must be positive');

assert(isequal(abs(b0), b0), 'b must be real-valued and non-negative');

if ~isnumeric(A) & isnumeric(At)
    error('If A is a function handle, then At must also be a function handle');
end

if ~isempty(opts.customx0)
    assert(isequal(size(opts.customx0), [n 1]), ...
        'customx0 must be a column vector of length n');
end

end

% Check that A and At are indeed ajoints of one another
function checkAdjoint(A, At, b)
y = randn(size(b));
Aty = At(y);
x = randn(size(Aty));
Ax = A(x);

innerProduct1 = Ax(:)'*y(:);
innerProduct2 = x(:)'*Aty(:);
error = abs(innerProduct1-innerProduct2)/abs(innerProduct1);
assert(error<1e-3 , ['Invalid measurement operator:  At is not the adjoint of A.  Error = ',num2str(error)]);
end

