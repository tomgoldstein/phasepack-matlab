%                     experimentTransMatrixWithSynthData.m
%
% Load an emprical tranmission matrix, and use it to acquire measurements
% from a synthetic random signal.
%
% Inputs:
%       n: The length of the signal
%       m: The number of measurements to sample from the transmission
%              matrix.  The m 'best' rows of the matrix will be used.
%       A_cached:  If the measurement matric 'A' has already been loaded,
%               then it can be recycled instead of re-loading.
%
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

%% -----------------------------START----------------------------------


function [A, b, xt, plotter] = experimentTransMatrixWithSynthData(n, m, A_cached)
    
    % Create some strings for dealing with different filenames according to the dataset chosen
    switch n
        case {16*16}
            measurementType = 'AmpSLM_16x16';
        case {40*40}
            measurementType = 'PhaseSLM_40x40';
        case {64*64}
            measurementType = 'AmpSLM_64x64';
        otherwise
            error('illegal dataSize: %d. It should be chosen from {256, 1600, 4096}',dataSize);
    end

    
    %% Load the transmission matrix, ground truth image, and measurements
    dataRoot = getFolderPath('data');       % The location of the 'data' folder
    % Load the data
    if numel(A_cached)>0 && isequal(size(A_cached),[m,n])
        A = A_cached;
    else
        try
            %disp(['    Loading transmission matrix (A_GS.mat), this may take several minutes']);
            load(strcat(dataRoot,'TransmissionMatrices/Coherent_Data/',measurementType,'/A_GS.mat')); 
        catch
            error(['You are missing the transmission matrix dataset.  To download it, Go to',...
                'https://rice.app.box.com/s/0c7thl2ck06mhaov1y3wnbzgjalf6ssq/folder/26509131173,'...
                'and download the "TransmissionMatrices" folder. Unzip the downloaded folder and',...
                'place it in the in "data" folder. See the user guide for details.']); 
        end
        % Only use the most accurate rows of the transmission matrix.  This is
        % determined by sorting the residuals computed during TM
        % reconstruction, and taking the rows with smallest residuals.
        [sortResid, ind] = sort(residual_vector,'ascend');
        assert(m<=size(A,1),['When using this dataset, you must choose m<',num2str(size(A,1))])
        ind = ind(1:m);
        A = A(ind,:);
    end
    
    % create synthetic measurements
    xt = randn(n,1)+1i*randn(n,1);
    b = abs(A*xt);
    
    % create a function to plot reconstructions
    plotter = @(x) plot(abs(x),abs(xt));
    
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        Utility functions                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Find the path for the specified folder by searching the current and
% parent directories
function path = getFolderPath(fname)
    d = dir;
    % check if the folder fname is in current path
    if any(strcmp(fname, {d([d.isdir]).name}))
        path = [fname,'/'];
    else
        % Look for folder fname in parent directory
        d = dir('../');
        if any(strcmp(fname, {d([d.isdir]).name}))
            path = ['../',fname,'/'];
        else
            error(['Unable to find path of folder: ',fname, ...
                '.  Make sure your current directory is the PhasePack root.']);
        end
    end
end


function plot(x,xt)
    scatter(abs(xt),abs(x));
    xlabel('true signal','fontsize',14);
    ylabel('recovered signal','fontsize',14);
end

