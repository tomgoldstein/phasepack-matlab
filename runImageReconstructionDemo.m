%%                   runImageReconstructionDemo.m
%
% This script will create phaseless measurements from a test image, and 
% then recover the image using phase retrieval methods.  We now describe 
% the details of the simple recovery problem that this script implements.
% 
%                         Recovery Problem
% This script loads a test image, and converts it to grayscale.
% Measurements of the image are then obtained by applying a linear operator
% to the image, and computing the magnitude (i.e., removing the phase) of 
% the linear measurements.
%
%                       Measurement Operator
% Measurement are obtained using a linear operator, called 'A', that 
% obtains masked Fourier measurements from an image.  Measurements are 
% created by multiplying the image (coordinate-wise) by a 'mask,' and then
% computing the Fourier transform of the product.  There are 8 masks,
% each of which is an array of binary (+1/-1) variables. The output of
% the linear measurement operator contains the Fourier modes produced by 
% all 8 masks.  The measurement operator, 'A', is defined as a separate 
% function near the end of the file.  The adjoint/transpose of the
% measurement operator is also defined, and is called 'At'.
%
%                         Data structures
% PhasePack assumes that unknowns take the form of vectors (rather than 2d
% images), and so we will represent our unknowns and measurements as a 
% 1D vector rather than 2D images.
%
%                      The Recovery Algorithm
% The image is recovered by calling the method 'solvePhaseRetrieval', and
% handing the measurement operator and linear measurements in as arguments.
% A struct containing options is also handed to 'solvePhaseRetrieval'.
% The entries in this struct specify which recovery algorithm is used.
%
% For more details, see the Phasepack user guide.
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

function runImageReconstructionDemo()

%% Specify the target image and number of measurements/masks
image = imread('data/logo.jpg'); % Load the image from the 'data' folder.
image = double(rgb2gray(image)); % convert image to grayscale
num_fourier_masks = 8;           % Select the number of Fourier masks


%% Create 'num_fourier_masks' random binary masks. Store them in a 3d array.
[numrows,numcols] = size(image); % Record image dimensions
random_vars = rand(num_fourier_masks,numrows,numcols); % Start with random variables
masks = (random_vars<.5)*2 - 1;  % Convert random variables into binary (+1/-1) variables

%% Compute phaseless measurements
% Note, the measurement operator 'A', and it's adjoint 'At', are defined
% below as separate functions
x = image(:);  % Convert the signal/image into a vector so PhasePack can handle it
b = abs(A(x));  % Use the measurement operator 'A', defined below, to obtain phaseless measurements.

%% Set options for PhasePack - this is where we choose the recovery algorithm
opts = struct;             % Create an empty struct to store options
opts.algorithm = 'twf';    % Use the truncated Wirtinger flow method to solve the retrieval problem.  Try changing this to 'Fienup'.
opts.initMethod = 'TruncatedSpectral';  % Use a truncated spectral method to generate an initial starting point for the solver  
opts.tol = 1e-3;           % The tolerance - make this smaller for more accurate solutions, or larger for faster runtimes
opts.verbose = 2;          % Print out lots of information as the solver runs (set this to 1 or 0 for less output)

%% Run the Phase retrieval Algorithm
fprintf('Running %s algorithm\n',opts.algorithm);
% Call the solver using the measurement operator 'A', its adjoint 'At', the
% measurements 'b', the length of the signal to be recovered, and the
% options.  Note, this method can accept either function handles or
% matrices as measurement operators.   Here, we use function handles
% because we rely on the FFT to do things fast.
[x, outs, opts] = solvePhaseRetrieval(@A, @At, b, numel(x), opts);

% Convert the vector output back into a 2D image
recovered_image = real(reshape(x,numrows,numcols));

% Phase retrieval can only recover images up to a phase ambiguity. 
% Let's apply a phase rotation to align the recovered image with the 
% original so it looks nice when we plot it.
rotation = sign(recovered_image(:)'*image(:));
recovered_image = recovered_image*rotation;

% Print some useful info to the console
fprintf('Image recovery required %d iterations (%f secs)\n',outs.iterationCount, outs.solveTimes(end));


%% Plot results
figure;
% Plot the original image
subplot(1,3,1);
imagesc(image);
title('Original Image');
% Plot the recovered image
subplot(1,3,2);
imagesc(real(recovered_image));
title('Recovered Image');
% Plot a convergence curve
subplot(1,3,3);
convergedCurve = semilogy(outs.solveTimes, outs.residuals);
set(convergedCurve, 'linewidth',1.75);
grid on;
xlabel('Time (sec)');
ylabel('Error');
title('Convergence Curve');
set(gcf,'units','points','position',[0,0,1200,300]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%               Measurement Operator Defined Below                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a measurement operator that maps a vector of pixels into Fourier
% measurements using the random binary masks defined above.
function measurements = A(pixels)
    % The reconstruction method stores iterates as vectors, so we need to
    % accept a vector input, and convert it back into an image.
    im = reshape(pixels,[numrows,numcols]);
    % Allocate space for all the measurements
    measurements= zeros(num_fourier_masks,numrows,numcols);
    % Loop over each mask, and obtain the Fourier measurements
    for m = 1:num_fourier_masks
        this_mask = squeeze(masks(m,:,:));
        measurements(m,:,:) = fft2(im.*this_mask);
    end
    % Convert results into vector format
    measurements = measurements(:);
end

% The adjoint/transpose of the measurement operator
function pixels = At(measurements)
    % The reconstruction method stores measurements as vectors, so we need 
    % to accept a vector input, and convert it back into a 3D array of 
    % Fourier measurements.
    measurements = reshape(measurements,[num_fourier_masks,numrows,numcols]);
    % Allocate space for the returned value
    im = zeros(numrows,numcols);
    for m = 1:num_fourier_masks
        this_mask = squeeze(masks(m,:,:));
        this_measurements = squeeze(measurements(m,:,:));
        im = im + this_mask.*ifft2(this_measurements)*numrows*numcols;
    end
    % Vectorize the results before handing them back to the reconstruction
    % method
    pixels = im(:);
end

end