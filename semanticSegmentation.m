function [L, idx, centroidsFinal, classes, parameters, superPixels, pixelsOwn, pixelsAdj, rgbImage, originalRgbImage] = semanticSegmentation(imagePath, m, myFilter2, classes, centroidsFinal, parameters, scaleSize)

%% Inicia o código e define pastas e imagens
% initCode;

%% Define flags
printResults = 1;
plotsIM = 0; % 0 DON'T PLOT IMAGES | 1 PLOT COLORED SUPERPIXELS
plotsCompare = 0; % 0 DON'T PLOT IMAGES | 1 PLOT FCM SEGMENTATION
myFilter = myFilter2; % 0 NO FILTER | 1 BILATERAL | 2 KUWAHARA | 3 ANISIOTROPIC
myColorSpace = 0; % 0 RGB | 1 HSV | 2 Lab | 3 sRGB
myFeatureExtractor = 4; % 1 COLOR | 2 LBP | 3 TEXTURE (LM FILTERS) | 4 COLOR+TEXTURE | 7 GLCM | 5 AUTO-ENCODER
myClassifier = 0; % 1 MSE | 2 SVM | 3 CANFIS | 4 COSINE | 5 ALL (COMPARISON)
myClusters = 1; % 0 MY CLASSES | 1 MY FCM | 2 MATLAB ANFIS
classifyMethod = 1; % 0 SUPERPIXELS | 1 CENTROIDS
useCropped = 0; % 0 UAV IMAGES | 1 CROPPED IMAGES
myPreprocess = 2; % 0 NO PREPROCESS | 1 HISTOGRAM MATCHING | 2 HISTOGRAM EQUALIZATION+MATCHING
usePCA = 1; % 0 ORIGINAL FEATURES | 1 SINGLE PCA SPACE | 2 INDIVIDUAL PCA SPACE
bestK = 0; % WRONG, FIX LATER - 0 USE PREDEFINED CLUSTER NUMBER | 1 SEARCH FOR BEST CLUSTER NUMBER
myNormalization = 3; % 0 NO NORMALIZATION | 1 FEATURE | 2 SUPERPIXEL | 3 FEATURE+SUPERPIXEL
combine = 0; % 0 ORIGINAL CLASSES | 1 COMBINED CLASSES
increaseSp = 1; % 0 ORIGINAL SP NUMBER | 1 INCREASED SP NUMBER
accCombinedAll = [];
myAccs = [];

%% Leitura e pre-processamento das imagens
readImagePath;
preProcessImages;
originalRgbImage = rgbImage;

%% Aplicando filtro
applyFilter;

%% Definindo espaço de cor
colorSpace;

%% Segmentando e achando K para superpixels
[x, y, ~] = size(rgbImage);
% K = 2*round((x/100) * (y/100));
K = round((x/100) * (y/100));
if increaseSp
    K = round((x/50) * (y/50));
end

% Segmentando super pixels
superpixelSegmentation;

%% Select super pixels
loadColors;
if myClusters == 0
    if (m == 1)
        getClasses;
    else
        getClassesValidate;
    end
end

%% Extraindo caracteristicas
extractFeatures;
myPCA;

%% Calculando K otimo

if m == 1
    calculaK;

    % Temporariamente deixando um valor fixo para agilizar experimentos
    if myClusters == 0
        KF = length(myClasses);
    end
    pixelsOri = pixels;
    centroidsFinal = zeros(KF, 1);

    % MATLAB ANFIS CLASSES
    if (myClusters == 2)
        maxIter = 100;
       inFis = genfis3(pixelsOri, classes, 'sugeno', length(unique(classesOri))); % gera sistema nebuloso
       opt = anfisOptions('InitialFIS', inFis, 'EpochNumber', maxIter, 'InitialStepSize', 0.02);
       [outFis, trainError] = anfis([pixelsOri classesOri], opt); % treina ANFIS
       for o = 1:length(outFis.Inputs)
            outFis.Inputs(o).Range = [0 1];
       end
    end

    % MY FCM CLASSES
    if (myClusters == 1)
        tic;

        Nexp = 1000;
        allU = cell(Nexp, 1);
        allJ = zeros(Nexp, 1);
        allCentroids = cell(Nexp, 1);
        allClasses = cell(Nexp, 1);
        allMetrics = zeros(Nexp, 3);
        for k = 1:Nexp
            [allClasses{k}, allU{k}, JcTemp, allCentroids{k}] = MyFuzzyMeans_opt(pixels, KF);
            allJ(k) = JcTemp(end);
            allMetrics(k, 1) = evalclusters(pixels, allClasses{k}, 'CalinskiHarabasz').CriterionValues; % maior, melhor
            allMetrics(k, 2) = evalclusters(pixels, allClasses{k}, 'DaviesBouldin').CriterionValues; % menor, melhor
            allMetrics(k, 3) = evalclusters(pixels, allClasses{k}, 'silhouette').CriterionValues; % maior, melhor
        end
        
        idxCal = find(allMetrics(:, 1) == max(allMetrics(:, 1)), 1);
        idxDav = find(allMetrics(:, 2) == min(allMetrics(:, 2)), 1);
        idxSil = find(allMetrics(:, 3) == max(allMetrics(:, 3)), 1);
        medianJc = abs(allJ - median(allJ));
        meanJc = abs(allJ - mean(allJ));
        idxJcMedian = find(medianJc == min(medianJc), 1);
        idxJcMean = find(meanJc == min(meanJc), 1);
%         idxJ = idxJcMedian;
%         idxJ = idxJcMean;
%         idxJ = find(allJ == min(allJ), 1);
        idxJ = idxCal;
        
        classes = allClasses{idxJ};
        U = allU{idxJ};
        J = allJ(idxJ);
        centroidsFinal = allCentroids{idxJ};
        
        clear allU allJ allCentroids allClasses
%         [classes, U, J, centroidsFinal] = MyFuzzyMeans_opt(pixelsOri, KF);
%         [~, idxU] = sort(U, 2);
%         classes = idxU(:, end);
        classesOri = double(U == max(U, [], 2));
%         centroidsReconstructed = ((parameters.pcaCoeffs(:,1:parameters.pcaN)*centroidsFinal')+parameters.pcaMean')';
%         textureCentroids = centroidsReconstructed(:, 85:end);
        fprintf('\nExecution time for fcm with %d classes: %f s\n', KF, toc);
    end

    fprintf('\nSaving cluster colors...\n');
    tic;
    saveOriginalColors;
    fprintf('Execution time for saving colors: %f s\n', toc);

    % MY TEST CLASSES
    if (myClusters == 0 || myClusters == 2)
        centroidsFinal = tempCentroids;
    end

    if plotsIM && myClassifier ~= 5
        colorSuperpixel;
        Lmask = zeros(size(L));
        for mx = 1:length(classes)
            Lmask(idx{mx}) = classes(mx);
        end
        figure; imagesc(Lmask);
    end

    % SVM 
    if (myClassifier == 2 || myClassifier == 5)
        fprintf('\nTraining multiclass SVM:');
        tic;
        model = svmtrain(classes, pixelsOri, '-q');
        parameters.model = model;
        fprintf('\nExecution time for training SVM: %f s\n', toc);
    end

    % ANFIS
    myMF;
    parameters.fis = inFis;
%         testAcc;
    if (myClassifier == 3 || myClassifier == 5)
        fprintf('\nTraining multiclass ANFIS:');
        tic;
        alfa = 0.01;
        nEp = 100;
        [~, cent, sig, pi, qi] = CANFIS(pixelsOri, classesOri, nEp, alfa, parameters.fis);
        parameters.c = cent;
        parameters.s = sig;
        parameters.p = pi;
        parameters.q = qi;
        fprintf('Execution time for multiclass (C)ANFIS: %f s\n', toc);
    end

    dbg = 1;

else
 %% Classificando imagens

 % Classificar superpixels
    if (classifyMethod == 0)
        if (myClassifier == 5)
            classifySuperpixelAll;
        else
            classifySuperpixel;
        end
    end

    % Clasificar centroides
    if (classifyMethod == 1)
        preprocessClusters;
        if (myClassifier == 5)
            classifyCentroidsAll;
        elseif (myClassifier == 0)
            classes = classesTemp;
        else
            classifyCentroids;
        end
    end
    if plotsIM && myClassifier ~= 5
        colorSuperpixel;
    end
end

fprintf('\n'); % printa um espaco para ficar mais bonito

end
 