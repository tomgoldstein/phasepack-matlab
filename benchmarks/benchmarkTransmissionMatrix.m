%                     benchmarkTransmissionMatrix.m
%
% Code adapted from Sudarshan Nagesh for the reconstruction of images from
% a transmission matrix.
%
% This script will reconstruct images from phaseless measurements that were
% acquired when imaging through a diffusive media.  The quality of the
% reconstruction will be computed by comparing to a known solution, and
% different algorithms are compared.
%   Three datasets are available with resolutions 16x16, 40x40, and 64x64.
% At each resolution, there are 5 sets of measurements, each corresponding
% to a different image.  The user must select which resolution and dataset 
% they want to run on.
%   As the script runs, all reconstructed images will be saved inside the
% 'benchmarkResults' folder.
%
% Note: The empirical dataset must be downloaded and installed into the
% 'data' directory for this to work. See the user guide.
%
%
%%I/O
%  Inputs:
%  imageSize            : Size of image to reconstruct.  Must be {16,40,64}
%  datasetSelection     : Choose which of the sets of measurements to use.
%                           Must be in {1,2,3,4,5}.
%  residualConstant     : Only use rows of the transmission matrix that
%                           had residuals less than this cutoff. Must be 
%                           between 0 and 1.  Recommend 0.4.
%  algorithms           : A cell-array of structs that describe algorithms.
%
%
%%References
%  Paper Title:   Coherent inverse scattering via transmission matrices:
%  Efficient phase retrieval algorithms and a public dataset.
%  Authors:       Christopher A. Metzler,  Manoj
%  K. Sharma,  Sudarshan Nagesh,  Richard G. Baraniuk,
%                 Oliver Cossairt,  Ashok Veeraraghavan
%
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

%% -----------------------------START----------------------------------


function benchmarkTransmissionMatrix(imageSize, datasetSelection, residualConstant, algorithms)
    
    %% Load the transmission matrix, and phaseless measurements
    [A,Xt,b] = loadData(residualConstant, datasetSelection, imageSize);
    
    % This string contains a text-based table with all the results, and is
    % printed after all calculations are done.
    % Results will also be printed incrementally as the benchmark runs.
    resultsTable = '\nAlgorithm  |  Measurement Error \n';
    
    %% Loop over each algorithm, and perform reconstruction and evaluation 
    for k=1:length(algorithms)
    	opts=algorithms{k};
        fprintf('Running Algorithm: %s\n',opts.algorithm);
        
        % Call a phase retrieval routine, and record the error of the result
        [reconError, measurementError] = reconstructSignal(A,Xt,b,imageSize,datasetSelection,opts);
        
        % Record the error in the table
        resultsTable = [resultsTable,...
                        pad(opts.algorithm,16),...
                        sprintf('%0.5f',measurementError),...
                        '\n'];
    end
    
    %% print the table
    fprintf(resultsTable);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                        Utility functions                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the measurement operator and measurements that correspond to the
% image size and dataset chosen by the user.  This method returns the
% measurement operator 'A' (a dense matrix), a ground truth solution
% 'Xt', and a set of measurements 'Y'.
function [A,Xt,b] = loadData(residualConstant, datasetSelection, imageSize)
    % Create some strings for dealing with different filenames according to the dataset chosen
    switch imageSize
        case {16}
            measurementType = 'AmpSLM_16x16';
        case {40}
            measurementType = 'PhaseSLM_40x40';
        case {64}
            measurementType = 'AmpSLM_64x64';
        otherwise
            error('illegal imageSize: %d. It should be chosen from {16, 40, 64}',imageSize);
    end
    
    % Make sure the dataset selection is valid
    if imageSize==16 || imageSize==40
        assert(datasetSelection <= 5);
    elseif imageSize==64
        assert(datasetSelection <= 6);
    end
    
    %% Load the transmission matrix, ground truth image, and measurements
    dataRoot = getFolderPath('data');       % The location of the 'data' folder
    try
        disp(['Loading transmission matrix (A_prVAMP.mat), this may take several minutes']);
        load(strcat(dataRoot,'TransmissionMatrices/Coherent_Data/',measurementType,'/A_prVAMP.mat')); 
        disp(['Loading test image (XH_test.mat)']);
        load(strcat(dataRoot,'/TransmissionMatrices/Coherent_Data/',measurementType,'/XH_test.mat'));
        disp(['Loading measurements (YH_squared_test.mat)']);
        load(strcat(dataRoot,'TransmissionMatrices/Coherent_Data/',measurementType,'/YH_squared_test.mat'));
    catch
        error(['You are missing the transmission matrix dataset.  To download it, Go to',...
            'https://rice.app.box.com/v/TransmissionMatrices/,'...
            'and download the "TransmissionMatrices" folder. Unzip the downloaded folder and',...
            'place it in the in "data" folder. See the user guide for details.']); 
    end
    % Only use the most accurate rows of the transmission matrix.  This is
    % determined by checking that the residuals found during row calcuation
    % are less than the specificed "residualConstant".
    findPixels = find(residual_vector<residualConstant);
    YH_squared_test = YH_squared_test(:,findPixels);
    
    A = A(findPixels,:);
    
    % Unpack the measurement data
    Y = double(YH_squared_test(datasetSelection,:));
    b = sqrt(Y);

    % Unpack the ground-truth image Xt
    % Note:  this image is real-valued
    Xt = double(XH_test(datasetSelection,:))';
    Xt = Xt/max(abs(Xt));
    % reshape original image
    originalImage = imresize(reshape(real(Xt),imageSize,imageSize),4);
    originalImageName = ['TM',num2str(imageSize),'-',num2str(datasetSelection),'-original'];
    imshow(originalImage,[]);  title(originalImageName,'fontsize',16); drawnow;
    % save the image to a folder
    saveBenchmarkImage(originalImage,['TM',num2str(imageSize)],originalImageName);
end

% Reconstruct a signal from a given measurement operator and set of
% measurements.  Compare the results to the ground-truth Xt, and report
% error.
function [reconError, measurementError] = reconstructSignal(A,Xt,Y,imageSize,datasetSelection,opts)
    fprintf('    Reconstructing image...')
    
    % Convenient variables
    n = imageSize^2;                    % Size of SLM
    imageSizeStr = num2str(imageSize);  % For TM16, this is '16'
    ind = datasetSelection;             % The index of the selected dataset
    
    % Solve the PR problem
    [X,outs,opts] = solvePhaseRetrieval(A,[],Y(:),n,opts);
   
    % Compute solution quality by solving a system to get the best affine
    % transform that maps the recovered solution onto the true solution
    Xrec = [X(:), ones(numel(X),1)];
    coeffs = Xrec\Xt(:);
    Xrec = Xrec*coeffs;
    
    % Reshape images into a square, and save them.
    realIm = imresize(reshape(real(Xrec),imageSize,imageSize),4);
    recoveredImageNameAbs = strcat('TM',imageSizeStr,'-',num2str(ind),'-',opts.algorithm,...
                                num2str(outs.iterationCount),'-real');
    imshow(realIm,[]); title(recoveredImageNameAbs,'fontsize',16); drawnow;
    realIm = realIm-min(realIm(:));
    realIm = realIm/max(realIm(:));
    saveBenchmarkImage(realIm,['TM',imageSizeStr],recoveredImageNameAbs);
    

    reconError = norm(Xt(:)-Xrec(:))/norm(Xt(:));
    measurementError = norm(abs(A*X)-Y(:))/norm(Y(:));
    %fprintf('    Relative reconstruction error = %s\n',reconError);
    fprintf('    Relative measurement error = %s\n',measurementError);
end
    
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

% Save an image to a sub-folder inside 'benchmarkResults'
function saveBenchmarkImage(image,folder,fname)
    bmRoot = getFolderPath('benchmarkResults');
    if ~exist([bmRoot,folder],'dir')
         mkdir([bmRoot,folder]);
    end
    fullName = strcat(bmRoot,folder,'/',fname,'.png');
    fprintf('    Saving reconstructed image: %s\n',fullName);
    try
        imwrite(image,fullName);
    catch
        error(['Unable to save image to file: ', fullName]);
    end
end

