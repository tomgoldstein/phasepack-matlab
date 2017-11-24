%                               generateOutput.m
%
% Generate output struct according to the convergence info recorded.
% 
% Inputs:
%    opts(struct)                 : consists of options by default of
%                                   chosen by user.
%                                   For details see the User Guide.
%    iter(integer)                : The total iterations a solver runs.
%    solveTimes(vector)           : consists of elapsed time at each
%                                   iteration. 
%    measurementErrors(vector)    : consists of measurement
%                                   error at each iteration. A single
%                                   measurement error at a certain
%                                   iteration is equal to
%                                   norm(abs(Ax)-b0)/norm(b0), where A is
%                                   the m x n measurement matrix or
%                                   function handle x is the n x 1
%                                   estimated signal at that iteration and
%                                   b0 is the m x 1 measurements.
%    reconErrors(vector)          : consists of relative reconstruction
%                                   error at each iteration.
%                                   A single reconstruction error at a
%                                   certain iteration is equal to
%                                   norm(xt-x)/norm(xt), where xt is the m
%                                   x 1 true signal, x is the n x 1
%                                   estimated signal at that iteration.
%    residuals(vector)            : consists of residuals at each
%                                   iteration.
%                                   Definition of a single residual depends
%                                   on the specific algorithm used see the
%                                   specific algorithm's file's header for
%                                   details.
% Outputs:
%    outs : A struct with convergence information
%    iterationCount(integer) : the number of
%                              iteration the algorithm runs. 
%    solveTimes(vector) : consists of elapsed (exist when
%                         recordTimes==true) time at each iteration.
%                         
%    measurementErrors(vector) : consists of the errors (exist when
%                                recordMeasurementErrors==true)   i.e.
%           norm(abs(A*x-b0))/norm(b0) at each iteration.
%           
%    reconErrors(vector): consists of the reconstruction (exist when
%                      recordReconErrors==true) errors i.e.
%                      norm(xt-x)/norm(xt) at each iteration.
%          
%    residuals(vector): consists of values that (exist when
%                       recordResiduals==true) will be compared with
%                       opts.tol for stopping condition checking.
%                       Definition varies across solvers.
%
%

% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

%% -----------------------------START----------------------------------


function outs = generateOutputs(opts, iter, solveTimes, measurementErrors, reconErrors, residuals)
    outs = struct;
    if ~isempty(solveTimes)
        outs.solveTimes = solveTimes(1:iter);
    end
    if ~isempty(measurementErrors)
        outs.measurementErrors = measurementErrors(1:iter);
    end
    if ~isempty(reconErrors)
        outs.reconErrors = reconErrors(1:iter);
    end
    if ~isempty(residuals)
        outs.residuals = residuals(1:iter);
    end
    outs.iterationCount = iter;
end