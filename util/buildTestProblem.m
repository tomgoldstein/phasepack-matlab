%                           buildTestProblem.m
%
% This function creates and outputs random generated data and measurements
% according to user's choice. It is invoked in test*.m in
% order to build a test problem.
%
% Inputs:
%   m(integer): number of measurements.
%   n(integer): length of the unknown signal.
%   isComplex(boolean, default=true): whether the signal and measurement
%     matrix is complex. isNonNegativeOnly(boolean, default=false): whether
%     the signal is real and non-negative.
%   dataType(string, default='gaussian'): it currently supports
%     ['gaussian', 'fourier'].
%
% Outputs:
%   A: m x n measurement matrix/function handle.
%   xt: n x 1 vector, true signal.
%   b0: m x 1 vector, measurements.
%   At: A n x m matrix/function handle that is the transpose of A.
%
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein
% Copyright (c) University of Maryland, 2017

function [A, xt, b0, At] = buildTestProblem(m, n, isComplex, isNonNegativeOnly, dataType)
if ~exist('isComplex','var') | isempty(isComplex)
    isComplex = true;
end
if ~exist('isNonNegativeOnly','var') | isempty(isNonNegativeOnly)
    isNonNegativeOnly = false;
end
if ~exist('dataType','var')
    dataType = 'Gaussian';
end

switch lower(dataType)
    case 'gaussian'
        A = (mvnrnd(zeros(1, n), eye(n)/2, m) + isComplex * 1i * ...
            mvnrnd(zeros(1, n), eye(n)/2, m));
        At = A';
        xt = (mvnrnd(zeros(1, n), eye(n)/2) + isComplex * 1i * ...
            mvnrnd(zeros(1, n), eye(n)/2))';
        b0 = abs(A*xt);
    case 'fourier'
        %  Define the Fourier measurement operator.
        %  The operator 'A' maps an n-vector into an m-vector, then
        %  computes the fft on that m-vector to produce m measurements.
        
        % rips first 'length' entries from a vector
        rip = @(x,length) x(1:length);
        A = @(x) fft([x;zeros(m-n,1)]);
        At = @(x) rip(m*ifft(x),n);     % transpose of FM
        xt = (mvnrnd(zeros(1, n), eye(n)/2) + isComplex * 1i * ...
            mvnrnd(zeros(1, n), eye(n)/2))';
        b0 = abs(A(xt)); % Compute the phaseless measurements
        
    otherwise
        error('invalid dataType: %s',dataType);
end
end