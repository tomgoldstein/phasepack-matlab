%                        experimentImage2D.m
%
% Create a Fourier measurement operator and measurements for reconstructing
% an image.
%
% Inputs:
% numMasks: number of random octanary Fourier masks. 
% imagePath: The path of the image to be recovered
%
% Outputs:
%  A     : A function handle: n1*n2 x 1 -> n1*n2*numMasks x 1. It returns 
%          A*x.
%  At    : The transpose of A
%  b0    : A n1*n2*numMasks x 1 real, non-negative vector consists of the 
%          measurements abs(A*x).
%  Xt    : The true signal - a vectorization of the image
%
%
%
% PhasePack by Rohan Chandra, Ziyuan Zhong, Justin Hontz, Val McCulloch,
% Christoph Studer, & Tom Goldstein 
% Copyright (c) University of Maryland, 2017

%% -----------------------------START----------------------------------


function [A, At, b0, Xt, plotter] = experimentImage2D(numMasks, imagePath)
    %%  Read Image 
    % Below Xt is n1 x n2 x 3; i.e. we have three n1 x n2 images, one for each of the 3 color channels  
    image = rgb2gray(imread([imagePath]));
    image = double(image);
    dims = size(image);
    L = numMasks;

    %% Make octanary masks and linear sampling operators  
    % Each mask has iid entries following the octanary pattern.
    % Entries have the form b = b1*b2 , where b1 is sampled
    % from {-1, 1, -i, i} with equal probability 1/4, and b2 from
    % { sqrt(1/2), sqrt(3)} with probability 4/5 and 1/5 respectively.
    b1 = [-1;1;-1i;1i];
    b2 = [repmat(sqrt(0.5),4,1); sqrt(3)];
    % Masks has size [n1,n2,L]
    masks = b1(randi(4,[dims,L])) .* b2(randi(5,[dims,L])); % Storage for L masks, each of dim n1 x n2
    
    %% Make linear operators that act on a vectorized image
    A = @(x) fourierMeasurementOperator(x, masks, dims);
    At = @(y) transposeOperator(y, masks, dims);
    
  
    Xt = image(:);
    b0 = abs(A(Xt(:)));
    
    % Set up plotting
    subplot(1,2,1);
    imagesc(image);
    title('original');
    subplot(1,2,2);
    
    plotter = @(x) imagesc(reshape(abs(x),dims));

end

function y = fourierMeasurementOperator(x, masks, dims)
    x = reshape(x, dims);   % image comes in as a vector.  Reshape to rectangle
    [n1,n2] = size(x);
    L = size(masks,3);              % get number of masks

    % Compute measurements
    copies = repmat(x,[1,1,L]);
    y = fft2(masks.*copies);
    y = y(:);
end


function x = transposeOperator(y, masks, dims)
    n1 = dims(1);
    n2 = dims(2);
    L = size(masks,3);              % get number of masks
    y = reshape(y, [n1,n2,L]);   % image comes in as a vector.  Reshape to rectangle
    
    x = n1*n2*ifft2(y).*conj(masks);
    x = sum(x,3);
    x = x(:);
end

