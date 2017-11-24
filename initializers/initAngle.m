%% ----------------------------initAngle.m-----------------------------

% Given the true solution of a phase retrieval problem, this method
% produces a random initialization that makes the specified angle with that
% solution.  This routine is meant to be used for benchmarking only; it
% enables the user to investigate the sensitivity of different methods on
% the accuracy of the initializer.
%  
%% I/O
%  Inputs:
%     Xtrue:  a vector
%     theta: an angle, specified in radians.
%% returns
%     x0: a random vector oriented with the specified angle relative to
%         Xtrue.
% 
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

%% -----------------------------START----------------------------------


function [x0] = initAngle(xt, theta)

% To get the correct angle, we add a perturbation to Xtrue.  
% Start by producing a random direction.
d = randn(size(xt));

% Orthogonalize the random direction
d = d - (d'*xt)/norm(xt)^2*xt;

% Re-scale to have same norm as the signal
d = d/norm(d)*norm(xt);

% Add just enough perturbation to get the correct angle
x0 = xt + d*tan(theta);

% Normalize 
x0 = x0/norm(x0)*norm(xt);

end

