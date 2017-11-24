%%                   runInitAngle.m

% Test that initAngle is able to produce an initialization vector with the
% correct angle with respect to the true solution.

%% ---------------------------------START-----------------------------


clc
clear
close all

% Parameters
n = 500;           % number of unknowns

Xt = randn(n,1);
X0 = initAngle(Xt,.5);

angle = acos( (X0'*Xt)/norm(X0)/norm(Xt)  );

fprintf('angle: %f\n', angle);