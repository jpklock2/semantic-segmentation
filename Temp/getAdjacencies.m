function [superPixels, pixelsOwn, pixelsAdj, L, classes, idx] = getAdjacencies(rgbImage, parameters, myFilter)

%% Inicia o código e define pastas e imagens
% initCode;

%% Define flags
printResults = 0;
plotsCompare = 0; % 0 DON'T PLOT IMAGES | 1 PLOT FCM SEGMENTATION
myFeatureExtractor = 4; % 1 COLOR | 2 LBP | 3 TEXTURE (LM FILTERS) | 4 COLOR+TEXTURE | 7 GLCM | 5 AUTO-ENCODER
myClusters = 1; % 0 MY CLASSES | 1 MY FCM | 2 MATLAB ANFIS
usePCA = 1; % 0 ORIGINAL FEATURES | 1 SINGLE PCA SPACE | 2 INDIVIDUAL PCA SPACE
myNormalization = 3; % 0 NO NORMALIZATION | 1 FEATURE | 2 SUPERPIXEL | 3 FEATURE+SUPERPIXEL
myClassifier = 2; % 1 MSE | 2 SVM | 3 CANFIS | 4 COSINE | 5 ALL (COMPARISON)
m = 2;

%% Aplicando filtro
originalRgbImage = rgbImage;
applyFilter;

%% Segmentando e achando K para superpixels
[x, y, ~] = size(rgbImage);
% K = 2*round((x/100) * (y/100));
% K = round((x/100) * (y/100));
K = round((x/50) * (y/50));

% Segmentando super pixels
superpixelSegmentation;

%% Extraindo caracteristicas
extractFeatures;
myPCA;

%% Classificando centroides
preprocessClusters;
classes = classesTemp;
% classifyCentroids;


if printResults
fprintf('\n'); % printa um espaco para ficar mais bonito
end


end

