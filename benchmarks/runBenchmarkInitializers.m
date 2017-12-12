%                   runBenchmarkInitializers.m
%
%   This function runs an assortment of initialization algorithms on test
%   data and plots the accuracy of each method as a function of the number
%   of samples used.
%
%   Note: for small values of m (the number of samples), the orthogonal
%   initializer might produce warnings because the spectral matrix
%   low-rank, and the smallest eigenvalue is not unique.
%


%% User-defined parameters.
 % The dimension of the signal to reconstruct. This must be 256, 1600, or 
 % 4096 when using the transmission matrix dataset.  It can be any positive
 % integer when using Gaussian data.
n = 256;       
% A list containing the numbers of samples for which each algorithm is run
m = n*[.5,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];  
numTrials = 10; % The number of random trials to run for each value of m

% Select the measurement matrix to use
% Valid options: {'transmissionMatrix','gaussian'}
%measurementOperator = 'transmissionMatrix';
measurementOperator = 'gaussian'; 

%% Allocate space for te results
results = zeros(numel(m),6);    % to store the reconstruction quality for each of the 6 algorithms
verbose = false;                % run all methods with text output turned off
A = [];                         % initialize the cached valeu of A for loading transmission matrices
%% Loop over the values of m, and store the performance of each method
fprintf('running trials...\n');
for m_index = 1:numel(m);
    fprintf('    m = %d\n',m(m_index));
    for trial = 1:numTrials
        
        % generate a random signal reconstruction problem
        switch lower(measurementOperator)
            case {'gaussian','synthetic'}
                isComplex = true;
                isNonNegativeOnly = false;
                [A, At, b0, xt, ~] = experimentGaussian1D(n, m(m_index), isComplex, isNonNegativeOnly);
            case {'tm','transmissionmatrix'}
                [A, b0, xt, plotter] = experimentTransMatrixWithSynthSignal(n, m(m_index), A);
                At = [];
            otherwise
                error(['Invalid dataset choice (',datatype,'): valid choices are "synthetic" and "transmissionMatrix"']);
        end
        
        % original spectral
        x = initSpectral(A,At,b0,n,false,true,verbose);
        results(m_index,1) = results(m_index,1)+abs(corr(x,xt));
          
        % truncated spectral
        x = initSpectral(A,At,b0,n,true,true,verbose);
        results(m_index,2) = results(m_index,2)+abs(corr(x,xt));
        
        % truncated amplitude flow
        x = initAmplitude(A,At,b0,n,verbose);
        results(m_index,3) = results(m_index,3)+abs(corr(x,xt));
        
        % weighted spectral
        x = initWeighted(A,At,b0,n,verbose);
        results(m_index,4) = results(m_index,4)+abs(corr(x,xt));
        
        % optimal spectral
        x = initOptimalSpectral(A,At,b0,n,true,verbose);
        results(m_index,5) = results(m_index,5)+abs(corr(x,xt));
        
        % the null initializer
        x = initOrthogonal(A, At, b0, n, verbose);
        results(m_index,6) = results(m_index,6)+abs(corr(x,xt));

    end

end
results = results/numTrials;
names = {'spectral', 'truncated', 'amplitude', 'weighted', 'optimal', 'orthogonal'}; 
autoplot(m,results/numTrials,names);
xlabel('number of samples','fontsize',16);
ylabel('cosine similarity','fontsize',16);
title(['initializer accuracy vs number of sample: n=',num2str(n)]);
l = legend('show', 'Location', 'northeastoutside');     % show the legend
set(l,'fontsize',16);


