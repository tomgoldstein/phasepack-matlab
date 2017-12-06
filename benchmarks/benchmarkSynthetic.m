%                           benchmarkSynthetic.m
%
% A general framework for benchmarking different phase retrieval algorithms
% using synthetic signals and either synthetic or real meaurement matrices.
%
% I/O
% Inputs
% xitem:  A string that describes the desired x-axis label of the plot that
%         that is produced by this benchmark.  Valid options are 
%         ['m/n', 'snr','iterations', 'time'].  When using a 2D image with
%         Fourier measurements, then 'masks' should be used instead of
%         'm/n' to conntrol the number of Fourier masks.
% xvalues: A list of scalar values for which the performance of each
%          algorithm is measured.
% yitem:  A string to appear as the y-axis label.  Value should be 
%         drawn from ['reconError', 'measurementError', 'correlation'];
% algorithms: a cell array of options structs,
%             where each struct is the same as the input parameter 'opts'
%             for solvePhaseRetrieval. See the example scripts for details.
% dataSet: The name of dataset used. Currently supported options are
%          ['1DGaussian', '2DImage', 'transmissionMatrix'].
%
% params: a struct of options containing the following fields:
%          verbose(boolean, default=false):      
%                   If true, the result of each trial will be reported.
%          numTrials(integer, default=1):        
%                   The number of trials each algorithm/dataset 
%                   combination will run.
%          policy(string, default='median'):     
%                   How to compute the final yvalue used for ploting from 
%                   the values one gets by running numTrials trials. It 
%                   currently supports
%                   ['median','mean','best','successRate'].
%          successConstant(real number,defualt=1e-5): 
%                   If the yvalue of the current  trial is less than this,
%                   the trial will be counted as a success. This parameter
%                   will only be used when policy='successRate'.
%          maxTime(positive real number,default=120) :  
%                   Max time allowed for a single algorithm.
%          recordSignals(boolean, default=false):
%                   Whether to record the recovered signal at each trial.
%
%
% Outputs
%          results : A 3D struct consisting of the errors(error
%                   metric is based on yitem chosen by the user) of all 
%                   trials of algorithm/dataset combinations. Coordinates 
%                   are (x-axis value, algorithm index, trial number).
% recoveredSignals: A 4D cell array consisting of the recovered signal at
%                   each trial for each algorithm. Coordinates are (x-axis
%                   value, algorithm index, current trial number, the index
%                   of the recovered signal).
%
% For more details, please look at the Phasepack User Guide.
%
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein
% Copyright (c) University of Maryland, 2017

%% -----------------------------START----------------------------------


% Main benchmark interface
function [finalResults, results, recoveredSignals] = benchmarkSynthetic(xitem, xvalues, yitem, algorithms, dataSet, params)

% Check if the inputs are valid. Note: the fields of algorithms will be checked by
% solvePhaseRetrieval implicitly.
checkValidityOfInputs(xitem, xvalues, yitem, dataSet);

% If params is not provided, create it.
if ~exist('params', 'var')
    params = struct;
end

% Provide default value/validate the params provided by the user.
params = manageOptionsForBenchmark(dataSet, params);

% Check if the params are valid for the dataSet chosen by the user.
params = checkValidityOfParams(dataSet, params, xitem);

% Get the labels for x, y axis, numTrials, policy, successConstant,
% and recordSignals.
% For details of what these represent, see the header in this file or the User Guide.
numTrials = params.numTrials;
policy = params.policy;
successConstant = params.successConstant;
recordSignals = params.recordSignals;

% create struct to store results of each trial.
results = zeros(length(xvalues), length(algorithms), numTrials);
if recordSignals
    recoveredSignals = cell(length(xvalues), length(algorithms), numTrials);
end

fprintf('Benchmarking on dataset %s with %s as x axis and the %s\n %s of %d trials as y axis...\n\n',dataSet,xitem,policy,yitem,numTrials);

% Loop over the xvalues
for p=1:length(xvalues)
    fprintf('Running trails: %s=%g\n',xitem,xvalues(p))
    % Loop over the algorithms
    for k=1:length(algorithms)
        [opts,params] = setupTrialParameters(algorithms{k},xitem,xvalues(p),dataSet,params);
        fprintf('  %s:',algorithms{k}.algorithm);
        if numTrials==1
            fprintf('\n');
        end
        % Loop over the random trials
        for q=1:numTrials
            if numTrials>1 && params.verbose==0
                fprintf('*');
            end
            
            % Create a random test problem with the right dimensions and SNR
            [A, At, b0, xt, plotter, params] = createProblemData(dataSet, params);
            n = numel(xt);
            opts.xt = xt;
            
            % Call the solver specified in opts
            startTime = tic;
            [x, outs, opts] = solvePhaseRetrieval(A, At, b0, n, opts);
            elapsedTime = toc(startTime);
            
            % Update plot
            plotter(x);
            title(sprintf('%s (%s=%s)',opts.algorithm,xitem,num2str(xvalues(p))),'fontsize',16);
            drawnow();

            % Calculate the value the user is looking for.
            yvalue = evalY(yitem, x, xt, A, b0);
            % Store the value for this trial in a table.
            results(p, k, q) = yvalue;
            if recordSignals
                recoveredSignals{p,k,q} = x;
            end

            % Report results of a trial if verbose is true.
            if params.verbose
                reportResult(opts.initMethod, opts.algorithm, xitem, xvalues(p),...
                    yitem, yvalue, elapsedTime, q);
            end
        end
        fprintf('\n');
    end
end

% Get final results in order to plot a comparison graph for algorithms
% at different x values.
finalResults = getFinalResults(results, policy, successConstant, yitem);

% Plot the comparison of performance graph among different chosen algorithms and
% initial methods combinations.
plotComparison(xitem,yitem,xvalues,algorithms,dataSet,finalResults,policy);

end



%% Helper functions

% Set the parameters needed to solve a phase retrieval problem with the
% specified dataset, xtime, and yitem.
function [opts,params] = setupTrialParameters(opts,xitem,xval,dataSet,params)
switch lower(xitem)
    case 'iterations'
        opts.maxIters = xval;    % Update algorithm's max iterations.
        opts.tol = 1e-10;              % Set algorithm's tolerance to a very
        % small number in order to make it
        % run maxIters iterations.
        opts.maxTime = params.maxTime; % Set a time limit for each single trial. 
        opts.isComplex = params.isComplex;
        opts.isNonNegativeOnly = params.isNonNegativeOnly;
    case 'm/n'
        if strcmp(lower(dataSet),'1dgaussian') || strcmp(lower(dataSet),'transmissionmatrix')
            params.m = round(params.n * xval); % Update m according m/n ratio.
        end
        opts.maxTime = params.maxTime; % Set a time limit for each single trial
    case 'snr'
        opts.maxTime = params.maxTime; % Set a time limit for each single trial.
        params.snr = xval;
    case 'time'
        opts.maxTime = xval;   % Set a time limit for each single trial.
        opts.tol = 1e-10;            % Set algorithm's tolerance to be very small
        opts.maxIters = 1e8;         % and algorithm's maxIters to be very large so it runs maxTime.
        opts.recordTimes = false;    % To save space, record nothing
        opts.recordResiduals = false;
    case 'masks'
        params.numMasks = xval;
        opts.maxTime = params.maxTime; % Set a time limit for each single trial
    case 'angle'
        opts.maxTime = params.maxTime; % Set a time limit for each single trial.
        opts.initMethod = 'angle';
        opts.initAngle = xval;
    otherwise
        % Raise an error if the given label for x axis is invalid.
        error(['invalid x label: ' xitem]);
end
end 


% Run a specified algorithm on a specific dataset and get results
% Inputs:
% dataSet and params are as defined in benchmarkPhaseRetrieval
% opts: a struct consists of options specified by user for running the algorithm
%       on a specified dataset.
% Outputs:
% x:  n x 1 vector, estimation of xt given by the PR algorithm the subrountine invokes.
% xt: n x 1 vector, the real unknown signal.
% A:  m x n matrix/function handle(n x 1 -> m x 1).
% b0: m x 1 vector, the measurement.
% opts: a struct consists of options finally used for running the algorithm on a specified
%       dataset.
function [A, At, b0, xt, plotter, params] = createProblemData(dataSet, params)
    switch lower(dataSet)
        case '1dgaussian'
            [A, At, b0, xt, plotter] = experimentGaussian1D(params.n, params.m,...
                params.isComplex, params.isNonNegativeOnly);
        case '2dimage'
            [A, At, b0, xt, plotter] = experimentImage2D(params.numMasks, params.imagePath);
        case 'transmissionmatrix'
            [A, b0, xt, plotter] = experimentTransMatrixWithSynthSignal(params.n, params.m, params.A_cached);
            params.A_cached = A;
            At = [];
        otherwise
            error('unknown dataset: %s\n',dataSet);
    end
    
    % Add noise to achieve specified SNR
    if params.snr ~= inf
        noise = randn(params.m,1); % create noise
        noise = noise/norm(noise)*norm(b0)/params.snr;
        b0 = max(b0+noise,0);
    end
    
end


% Calculate how good the solution is. Use the metric specified by yitem.
function yvalue=evalY(yitem, x, xt, A, b0)
    %  Compute optimal rotation
    switch lower(yitem)
        case 'reconerror'
            % solve for least-squares solution:  alpha*x = xt
            alpha = (x(:)'*xt(:))/(x(:)'*x(:));
            x = alpha*x;
            yvalue = norm(xt(:)-x(:))/norm(xt(:));
        case 'measurementerror'
            % Transform A into function handle if A is a matrix
            if isnumeric(A)
                At = @(x) A'*x;
                A = @(x) A*x;
            end
            yvalue = norm(abs(A(x))-b0(:))/norm(b0(:));
        case 'correlation'
            yvalue = abs(x'*xt/norm(x)/norm(xt));
        otherwise
            error(['invalid y label: ' yitem]);
    end
end

% Get final results by averaging across all trials. 
% The possible policies are the following:
% mean: take the mean of all trials.
% best: take the best of all trials. If yitem=='reconerror' or 'measurementerror',
%       min value will be taken; If yitem=='correlation', max value will be taken.
% median: take the median of all trials.
% successrate: A success rate will be calculated. If yitem=='reconerror' or
%              'measurementerror', it is the percentage of values that are
%              smaller than the successConstant. If yitem=='correlation', it is
%              the percentage of values that are larger than the successConstant.
% The input struct results has size length(xvalues) x length(algorithms) x numTrials
% The output struct finalResults has size length(xvalues) x length(algorithms)
function finalResults = getFinalResults(results, policy, successConstant, yitem)
    switch lower(policy)
        case 'mean'
            finalResults = mean(results,3);
        case 'best'
            switch lower(yitem)
                case {'reconerror', 'measurementerror'}
                    finalResults = min(results,[],3);
                case 'correlation'
                    finalResults = max(results,[],3);
                otherwise
                    error('invalid yitem: %s',yitem);
            end
        case 'median'
            finalResults = median(results,3);
        case 'successrate'
            switch lower(yitem)
                case {'reconerror', 'measurementerror'}
                    finalResults = mean(results<successConstant,3);
                case 'correlation'
                    finalResults = mean(results>successConstant,3);
                otherwise
                    error('invalid yitem: %s',yitem);
            end
        otherwise
            error('Invalid policy: %s', policy);
    end
end

% Plot a performance curve for each algorithm
function plotComparison(xitem,yitem,xvalues,algorithms,dataSet,finalResults,policy)
    algNames = {};
    % get the labels to appear in the legend
    for k=1:length(algorithms)
        if isfield(algorithms{k},'label') && numel(algorithms{k}.label)>0             % use user-specified label if available
            algNames{k} =  algorithms{k}.label;
        else
            algNames{k} = algorithms{k}.algorithm;  % otherwise use the algorithm name
        end
    end

    autoplot(xvalues,finalResults,algNames);

    title(char(strcat({xitem},{' VS '},{yitem},{' on '},{dataSet})));                         % title for plot
    xlabel(xitem);                                       % x-axis label
    ylabel(char(strcat({policy},{'('},{yitem},{')'})));  % y-axis label
    ax = gca;                                           % create an axes object
    ax.FontSize = 12;                                   % adjust font size on the axes
    legend('show', 'Location', 'northeastoutside');     % show the legend
end


% Report the result at the end of each trial if verbose is true.
function reportResult(initMethod, algorithm, xitem, xvalue, yitem, yvalue, time, currentTrialNum)
    if currentTrialNum==1
        fprintf('\n');
    end
    fprintf('Trial: %d, initial method: %s, algorithm: %s, %s: %d, %s: %f, time: %f\n',...
        currentTrialNum, initMethod, algorithm, xitem, xvalue, yitem, yvalue, time);
end


% Check if the inputs are valid.
function checkValidityOfInputs(xitem, xvalues, yitem, dataSet)
    assert(~isempty(xvalues), 'The list xvalues must not be empty.');

    % Check if yitem is a valid choice
    yitemList = {'reconerror', 'measurementerror', 'correlation'};
    checkIfInList('yitem', yitem, yitemList);

    % Check if dataSet chosen supports xitem
    xitem = lower(xitem);
    switch lower(dataSet)
        case '1dgaussian'
            supportedList = {'iterations','m/n','snr','time','angle'};
        case 'transmissionmatrix'
            supportedList = {'iterations','m/n','snr','time','angle'};
        case '2dimage'
            supportedList = {'masks','iterations','angle'};
        otherwise
            error('unknown dataset: %s\n',dataSet);
    end
    checkIfInList('xitem',xitem,supportedList,strcat('For ',lower(dataSet),', '));
end


% Check if the params are valid
function params = checkValidityOfParams(dataSet, params, xitem)
    switch lower(dataSet)
        case '1dgaussian'
            assert(isfield(params,'n'),'User must specify signal dimension in params.n');
            checkIfNumber('params.n',params.n);
            if ~strcmp(xitem,'m/n')
                assert(isfield(params,'m'), ['User must specify the number of measurements in params.m when using xlabel ',xitem]);
            end
        case '2dimage'
            assert(~strcmp(xitem,'m/n'),'xlabel m/n is not supported when using 2D images.  Instead use "masks"');
            if ~(xitem=='masks')
                assert(exist('params.numMasks'), ['User must specify the number of Fourier masks in params.numMasks when using xlabel ',xitem]);
            end
        case 'transmissionmatrix'
            params.A_cached = [];
            assert(isfield(params,'n'),'User must specify signal dimension in params.n.  Options are {256,1600,4096} for the tranmissionMatrix dataset.');
            assert(params.n==256 || params.n==1600 || params.n==4096,...
                ['Invalid choice (',num2str(params.n),') ofr params.n for the transmissionMatrix dataset. Valid choices are {256,1600,4096}.'])
        otherwise
            error('unknown dataset: %s\n',dataSet);
    end

    if ~isfield(params,'snr')
        params.snr = Inf;
    end

end
