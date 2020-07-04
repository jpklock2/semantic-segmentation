%% Inicia o c�digo e define pastas e imagens
initCode;

%% Define flags
plotsIM = 1; % 0 DON'T PLOT IMAGES | 1 PLOT COLORED SUPERPIXELS
plotsCompare = 0; % 0 DON'T PLOT IMAGES | 1 PLOT FCM SEGMENTATION
myFilter = 1; % 0 NO FILTER | 1 BILATERAL | 2 KUWAHARA | 3 ANISIOTROPIC
myColorSpace = 0; % 0 RGB | 1 HSV | 2 Lab | 3 sRGB
myFeatureExtractor = 4; % 1 COLOR | 2 LBP | 3 TEXTURE (LM FILTERS) | 4 COLOR+TEXTURE | 7 GLCM | 5 AUTO-ENCODER
myClassifier = 2; % 1 MSE | 2 SVM | 3 CANFIS | 4 COSINE | 5 ALL (COMPARISON)
myClusters = 0; % 0 MY CLASSES | 1 MY FCM | 2 MATLAB ANFIS
classifyMethod = 1; % 0 SUPERPIXELS | 1 CENTROIDS
useCropped = 1; % 0 UAV IMAGES | 1 CROPPED IMAGES
myPreprocess = 0; % 0 NO PREPROCESS | 1 HISTOGRAM MATCHING | 2 HISTOGRAM EQUALIZATION+MATCHING
usePCA = 1; % 0 ORIGINAL FEATURES | 1 SINGLE PCA SPACE | 2 INDIVIDUAL PCA SPACE
bestK = 0; % 0 USE PREDEFINED CLUSTER NUMBER | 1 SEARCH FOR BEST CLUSTER NUMBER
myNormalization = 3; % 0 NO NORMALIZATION | 1 FEATURE | 2 SUPERPIXEL | 3 FEATURE+SUPERPIXEL
combine = 0; % 0 ORIGINAL CLASSES | 1 COMBINED CLASSES
increaseSp = 0; % 0 ORIGINAL SP NUMBER | 1 INCREASED SP NUMBER
accCombinedAll = [];
myAccs = [];

 for m = 1:size(imageNames, 1) % iterando pelas imagens
     
     %% Leitura e pre-processamento das imagens
     
%      if (m == 1 || m == 2 || m == 3 || m == 8 || m == 18 || m == 31 || m == 47 || m == 58 || m == 79 || m == 206 || m == 294)
    
     pause(0.1);
     readImage;
     preProcessImages;
%      cropImages;
       
     %% Aplicando filtro
     applyFilter;
     
     %% Definindo espa�o de cor
     colorSpace;

    %% Segmentando e achando K para superpixels
    [x, y] = size(rgbImage);
    K = round((x/100) * (y/100));
    if increaseSp
        K = round((x/50) * (y/50));
    end
    
    % Segmentando super pixels
    superpixelSegmentation;
    
    %% Select super pixels
    loadColors;
    if (m == 1)
        getClasses;
    else
%         selectSuperpixels;
        getClassesValidate;
    end
%     combineClasses;
     
    %% Extraindo caracteristicas
    extractFeatures;
    myPCA;
    
    %% Calculando K otimo
    
    if m == 1
        calculaK;
    
        % Temporariamente deixando um valor fixo para agilizar experimentos
        if (myClusters == 0) 
            KF = length(myClasses);
        end
        pixelsOri = pixels;
        centroidsFinal = zeros(KF, 1);
%         KF = 4;
        
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
            [~, U, J, centroidsFinal] = MyFuzzyMeans_opt(pixelsOri, KF);
            [~, idxU] = sort(U, 2);
            classes = idxU(:, end);
            classesOri = double(U == max(U, [], 2));
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
        end
        
        % SVM 
        if (myClassifier == 2 || myClassifier == 5)
            fprintf('\nTraining multiclass SVM:');
            tic;
%             classesOriTemp = classesOri;
%             classesOriTemp(classesOriTemp == 0) = -1;
            model = svmtrain(classes, pixelsOri, '-q');
            fprintf('\nExecution time for training SVM: %f s\n', toc);
        end
        
        % ANFIS
        myMF;
%         testAcc;
        if (myClassifier == 3 || myClassifier == 5)
            fprintf('\nTraining multiclass ANFIS:');
            tic;
            alfa = 0.01;
            nEp = 100;
%             [~, cent, sig, pi, qi] = MANFIS(pixelsOri, classesOri, nEp, alfa, inFis);
            [~, cent, sig, pi, qi] = CANFIS(pixelsOri, classesOri, nEp, alfa, inFis);
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
            else
                classifyCentroids;
            end
        end
        if plotsIM && myClassifier ~= 5
            colorSuperpixel;
        end
    end

%     end
    
 end
% save(['acc_F' num2str(myFilter) '_CS' num2str(myColorSpace) '_FE' num2str(myFeatureExtractor) ...
%       '_C' num2str(myClassifier) '_A' num2str(useAnfis) '_CM' num2str(classifyMethod) '.mat'], 'myAccs');
% fprintf('\n Mean Accuracy = %.2f%%\n', mean(myAccs)); 
fprintf('\n'); % printa um espaco para ficar mais bonito
 