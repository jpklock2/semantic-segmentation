% clear all;
% clc;
function [edge_final] = Canny_v0(img, filter, kValue, expe, minIdx)
%Input image
% img = imread ('House.jpg');
%Show input image
% figure, imshow(img);
% img = rgb2gray(img);
% img = double (img);

%% MATLAB imp

sigma = sqrt(2);

% Magic numbers
PercentOfPixelsNotEdges = .7; % Used for selecting thresholds
ThresholdRatio = .4;          % Low thresh is this fraction of the high.

% Determine filter length
filterExtent = ceil(4*sigma);
x = -filterExtent:filterExtent;

% Create 1-D Gaussian Kernel
c = 1/(sqrt(2*pi)*sigma);
gaussKernel = c * exp(-(x.^2)/(2*sigma^2));

% Normalize to ensure kernel sums to one
gaussKernel = gaussKernel/sum(gaussKernel);

% Create 1-D Derivative of Gaussian Kernel
derivGaussKernel = gradient(gaussKernel);

% Normalize to ensure kernel sums to zero
negVals = derivGaussKernel < 0;
posVals = derivGaussKernel > 0;
derivGaussKernel(posVals) = derivGaussKernel(posVals)/sum(derivGaussKernel(posVals));
derivGaussKernel(negVals) = derivGaussKernel(negVals)/abs(sum(derivGaussKernel(negVals)));

% Compute smoothed numerical gradient of image I along x (horizontal)
% direction. GX corresponds to dG/dx, where G is the Gaussian Smoothed
% version of image I.
GX = imfilter(img, gaussKernel', 'conv', 'replicate');
GX = imfilter(GX, derivGaussKernel, 'conv', 'replicate');

% Compute smoothed numerical gradient of image I along y (vertical)
% direction. GY corresponds to dG/dy, where G is the Gaussian Smoothed
% version of image I.
GY = imfilter(img, gaussKernel, 'conv', 'replicate');
GY  = imfilter(GY, derivGaussKernel', 'conv', 'replicate');

% Calculate Magnitude of Gradient
magGrad = hypot(GX, GY);

% Normalize for threshold selection
magmax = max(magGrad(:));
if magmax > 0
    magGrad = magGrad / magmax;
end

[m,n] = size(magGrad);

% Select the thresholds
counts=imhist(magGrad, 64);
highThresh = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n,...
    1,'first') / 64;
lowThresh = ThresholdRatio*highThresh;

%% Other imp

%Value for Thresholding
% T_Low = 0.075;
% T_High = 0.175;
T_Low = highThresh;
T_High = lowThresh;

%Convolution of image by Gaussian Coefficient
if filter == 1
%Gaussian Filter Coefficient
B = [2, 4, 5, 4, 2; 4, 9, 12, 9, 4;5, 12, 15, 12, 5;4, 9, 12, 9, 4;2, 4, 5, 4, 2 ];
B = 1/159.* B;
A = conv2(img, B, 'same');

elseif filter == 2
k = 3;
A = Kuwahara(img,4*k+1);

elseif filter == 3
% patch = imcrop(img,[34,71,60,55]);
patch = minIdx;
patchVar = std2(img(patch)).^2;
% disp(patchVar);
DoS2 = 4*patchVar;
A = imbilatfilt(img,DoS2,4);

elseif filter == 4
lambda = 5;
iter = 5;
p = 0.2;
eps = 0.0001;
A = ILS_LNorm(img, lambda, p, eps, iter);

elseif filter == 5
[gradThresh,numIter] = imdiffuseest(img,'ConductionMethod','exponential');
diffImg = imdiffusefilt(img,'ConductionMethod','exponential', ...
'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
A = diffImg;

end

% figure; montage({img, A, A2, A3});

%Filter for horizontal and vertical direction
KGx = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
KGy = [1, 2, 1; 0, 0, 0; -1, -2, -1];

%Convolution by image by horizontal and vertical filter
if expe == 0 || expe == 1
Filtered_X = conv2(A, KGx, 'same');
Filtered_Y = conv2(A, KGy, 'same');
else
Filtered_X = imfilter(A, KGx', 'conv', 'replicate');
Filtered_Y = imfilter(A, KGy', 'conv', 'replicate');
end

%Calculate directions/orientations
arah = atan2 (Filtered_Y, Filtered_X);
arah = arah*180/pi;

pan=size(A,1);
leb=size(A,2);

%Adjustment for negative directions, making all directions positive
for i=1:pan
    for j=1:leb
        if (arah(i,j)<0) 
            arah(i,j)=360+arah(i,j);
        end
    end
end

arah2=zeros(pan, leb);

%Adjusting directions to nearest 0, 45, 90, or 135 degree
for i = 1  : pan
    for j = 1 : leb
        if ((arah(i, j) >= 0 ) && (arah(i, j) < 22.5) || (arah(i, j) >= 157.5) && (arah(i, j) < 202.5) || (arah(i, j) >= 337.5) && (arah(i, j) <= 360))
            arah2(i, j) = 0;
        elseif ((arah(i, j) >= 22.5) && (arah(i, j) < 67.5) || (arah(i, j) >= 202.5) && (arah(i, j) < 247.5))
            arah2(i, j) = 45;
        elseif ((arah(i, j) >= 67.5 && arah(i, j) < 112.5) || (arah(i, j) >= 247.5 && arah(i, j) < 292.5))
            arah2(i, j) = 90;
        elseif ((arah(i, j) >= 112.5 && arah(i, j) < 157.5) || (arah(i, j) >= 292.5 && arah(i, j) < 337.5))
            arah2(i, j) = 135;
        end
    end
end

% figure, imagesc(arah2); colorbar;

%Calculate magnitude
magnitude = (Filtered_X.^2) + (Filtered_Y.^2);
magnitude2 = sqrt(magnitude);

% Normalize for threshold selection
magmax = max(magnitude2(:));
if magmax > 0
    magGrad = magGrad / magmax;
end

[m,n] = size(magGrad);

% Select the thresholds
counts=imhist(magGrad, 64);
highThresh2 = find(cumsum(counts) > PercentOfPixelsNotEdges*m*n,...
    1,'first') / 64;
lowThresh2 = ThresholdRatio*highThresh2;

BW = zeros (pan, leb);

%Non-Maximum Supression
for i=2:pan-1
    for j=2:leb-1
        if (arah2(i,j)==0)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i,j+1), magnitude2(i,j-1)]));
        elseif (arah2(i,j)==45)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j-1), magnitude2(i-1,j+1)]));
        elseif (arah2(i,j)==90)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j), magnitude2(i-1,j)]));
        elseif (arah2(i,j)==135)
            BW(i,j) = (magnitude2(i,j) == max([magnitude2(i,j), magnitude2(i+1,j+1), magnitude2(i-1,j-1)]));
        end
    end
end

BW = BW.*magnitude2;
% figure, imshow(BW);

if expe == 1 || expe == 3
T_Low = highThresh2;
T_High = lowThresh2;
end

%Hysteresis Thresholding
T_Low = T_Low * max(max(BW)) * kValue;
T_High = T_High * max(max(BW)) * kValue;

T_res = zeros (pan, leb);

for i = 1  : pan
    for j = 1 : leb
        if (BW(i, j) < T_Low)
            T_res(i, j) = 0;
        elseif (BW(i, j) > T_High)
            T_res(i, j) = 1;
        %Using 8-connected components
        elseif ( BW(i+1,j)>T_High || BW(i-1,j)>T_High || BW(i,j+1)>T_High || BW(i,j-1)>T_High || BW(i-1, j-1)>T_High || BW(i-1, j+1)>T_High || BW(i+1, j+1)>T_High || BW(i+1, j-1)>T_High);
            T_res(i,j) = 1;
        end
    end
end

edge_final = logical(T_res);
% edge_final = uint8(T_res.*255);
%Show final edge detection result
% figure, imshow(edge_final);
end