%                                   experimentGaussian1D.m
%
% Produce a random signal reconstruction problem with a random Gaussian
% signal and a random Gaussian measurement matrix.
%
% Inputs:
%  n                : number of unknowns.
%  m                : number of measurements.
%  isComplex        : If the signal is complex.
%  isNonNegativeOnly: If the signal is non-negative.
%  SNR              : Signal-Noise ratio.
%
% Outputs:
%  A     : A function handle: nx1 -> mx1. It returns 
%          A*x.
%  At    : The transpose of A
%  b0    : A mx1 real, non-negative vector consists of the 
%          measurements abs(A*x).
%  xt    : The true signal.
%
%
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017


function [A, At, b0, xt, plotter] = experimentGaussian1D(n, m, isComplex, isNonNegativeOnly)
    
    % Build the test problem
    xt = randn(n, 1)+isComplex*randn(n, 1)*1i; % true solution

    if isNonNegativeOnly                       
        xt = abs(xt);
    end
    
    A = randn(m, n)+isComplex*randn(m, n)*1i; % matrix
    b0 = abs(A*xt);                           % measurements
    
    % Build function handles for the problem
    At = @(x) A'*x;
    A = @(x) A*x;
   
    
    plotter = @(x) plot(x,xt);

end

function plot(x,xt)
    scatter(abs(xt),abs(x));
    xlabel('true signal','fontsize',14);
    ylabel('recovered signal','fontsize',14);
end

