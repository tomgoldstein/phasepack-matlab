%                       plotErrorConvergence.m
%
% This function plots some convergence curve according to the values of
% options in opts specified by user. It is used in all the test*.m scripts.
% Specifically,
% If opts.recordReconErrors is true, it plots the convergence curve of
% reconstruction error versus the number of iterations. 
% If opts.recordResiduals is true, it plots the convergence curve of
% residuals versus the number of iterations. 
% The definition of residuals is algorithm specific. For details, see the
% specific algorithm's solve*.m file.
% If opts.recordMeasurementErrors is true, it plots the convergence curve
% of measurement errors. 
%
% Inputs are as defined in the header of solvePhaseRetrieval.m.
% See it for details.
%   
% 
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

function plotErrorConvergence(outs, opts)
    % Plot the error convergence curve
    if opts.recordReconErrors==true
        figure;
        semilogy(outs.reconErrors);
        xlabel('Iterations');
        ylabel('ReconErrors');
        title(strcat('Convergence curve:',{' '},opts.algorithm));
    end
    if opts.recordResiduals==true
        figure;
        semilogy(outs.residuals);
        xlabel('Iterations');
        ylabel('Residuals');
        title(strcat('Convergence curve:',{' '},opts.algorithm));
    end
    if opts.recordMeasurementErrors==true
        figure;
        semilogy(outs.measurementErrors);
        xlabel('Iterations');
        ylabel('MeasurementErros');
        title(strcat('Convergence curve:',{' '},opts.algorithm));
    end
end