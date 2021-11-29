function [modelSvm, modelTree] = createMasseliClassifier(rgbImage)

    plotFigs = 0;

    %% Preprocessing Image
    rgbImage = imresize(rgbImage, 0.25);

    labImage = rgb2lab(rgbImage);
    L = labImage(:,:,1)/100;
    L = adapthisteq(L);
    labImage(:,:,1) = L*100;
    rgbImage = lab2rgb(labImage);

    load('../CnnComparison/Images/Train/Classes/myClasses.mat');

    loadColors;
    myColorsClasses = ["Grass", "Forest", "Soil", "River"];
    myColors = myColors([2 4 8 1], :);
    
    % combining classes
    myClasses{1} = [myClasses{1}; myClasses{2}; myClasses{3}];
    myClasses{2} = [myClasses{4}; myClasses{5}; myClasses{12}];
    myClasses{3} = [myClasses{7}; myClasses{8}; myClasses{11}; myClasses{10}; myClasses{9}];
    myClasses{4} = [myClasses{6}; myClasses{13}; myClasses{14}];
    myClasses = myClasses(1:4);

    idx = label2idx(L2);
    Lmask = zeros(size(L2));
    for mx = 1:length(myClasses)
        Lmask(ismember(L2, myClasses{mx})) = mx;
    end
    
    %% Feature Extraction
    
    % pegando retangulo do meio
    x0 = 1;
    y0 = 1;
    x1 = size(rgbImage, 2);
    y1 = size(rgbImage, 1);

    usefulGray = rgb2gray(rgbImage);
    pad = 0;
    usefulGrayPad = padarray(usefulGray, [pad pad], 0, 'both');
    
    tY = y1-y0+2*pad;
    tX = x1-x0+2*pad;
    delta = 170;
    
    if plotFigs
        figure;
        imshow(usefulGray);
        hold on;
    end
    fts = [];
    labs = [];
    bestPoints = {[]};
    
    totalY = length(pad+1:delta:tY-delta-pad);
    totalX = length(pad+1:delta:tX-delta-pad);
    
    runFt = 1;
    if runFt
    
    cntY = 1;
    for ponY = pad+1:delta:tY-delta-pad
        fprintf('\nLinha %d/%d\n', cntY, totalY);
        cntX = 1;
        for ponX = pad+1:delta:tX-delta-pad
            fprintf('Coluna %d/%d\n', cntX, totalX);
            currImg = usefulGrayPad(ponY-pad:ponY-1+delta+pad, ponX-pad:ponX-1+delta+pad);
            currLabel = Lmask(ponY-pad:ponY-1+delta-pad, ponX-pad:ponX-1+delta-pad);
            
            [uni,~,ix] = unique(currLabel);
            acc = accumarray(ix,1);
            lab = uni(find(acc == max(acc), 1));
            
            Lmask(ponY-pad:ponY-1+delta-pad, ponX-pad:ponX-1+delta-pad) = lab;
            
%             points = detectORBFeatures(currImg);
            points = detectSURFFeatures(currImg);
            %points = detectORBFeatures(currImg, 'ScaleFactor', 1.5, 'NumLevels', 2);
            if isempty(points)
                fprintf('VAZIO!');
            end
            
            bestIdx = find(points.Metric == max(points.Metric), 1);
            bestP = points(bestIdx);
            bestP.Location(1) = bestP.Location(1) + ponX - 2*pad;
            bestP.Location(2) = bestP.Location(2) + ponY - 2*pad;
            
            tempFts = extractFeatures(currImg, points);
%             fts = [fts; tempFts.Features(bestIdx, :)];
            fts = [fts; tempFts(bestIdx, :)];
            bestPoints = [bestPoints; {bestP}];
            labs = [labs; lab];
%             figure
%             imshow(usefulGray)
%             hold on
            %plot(bestP);
%             hold off
            cntX = cntX+1;
        end
        cntY = cntY+1;
    end
    
    fts = double(fts);
    bestPoints(1) = [];
    save messeli3.mat bestPoints fts labs Lmask
    
    else
       load('messeli3.mat'); 
    end
    
    if plotFigs
        for i = 1:length(bestPoints)
            plot(bestPoints{i});
        end
    end
%     J = insertMarker(I,corners,'circle');
%     imshow(J)
    if plotFigs
        
        remainClasses = unique(Lmask);
        LmaskRed = zeros(size(L2));
        LmaskGren = zeros(size(L2));
        LmaskBlue = zeros(size(L2));
        LmaskColor = zeros(size(L2, 1), size(L2, 2), 3);
        LmaskCategorical = strings(size(L2));
        sizesClasses = zeros(1, length(myClasses));
        for mx = 1:length(remainClasses)
            LmaskRed(ismember(Lmask, remainClasses(mx))) = myColors(mx, 1);
            LmaskGren(ismember(Lmask, remainClasses(mx))) = myColors(mx, 2);
            LmaskBlue(ismember(Lmask, remainClasses(mx))) = myColors(mx, 3);
            LmaskCategorical(ismember(Lmask, remainClasses(mx))) = myColorsClasses(mx);
            sizesClasses(mx) = length(myClasses{mx});
        end
        LmaskColor(:, :, 1) = LmaskRed;
        LmaskColor(:, :, 2) = LmaskGren;
        LmaskColor(:, :, 3) = LmaskBlue;
        LmaskCategorical = categorical(LmaskCategorical, myColorsClasses);
    
    end
    
%     fts = fts./255;
    modelTree = TreeBagger(100,fts,labs,'Method','classification','OOBPrediction','On','NumPredictorsToSample',100);%,'NumPredictorsToSample',100);
%     modelTree = fitcensemble(fts,labs);
    classesT = str2num(cell2mat(predict(modelTree, fts)));
%     classesT = predict(modelTree, fts);
    accuracyTree = sum(classesT == labs)/length(labs);
    
    modelSvm = svmtrain(labs, fts, '-q');
    [classesS, accuracySvm, prob_estimates] = svmpredict(rand([size(fts, 1) 1]), fts, modelSvm, '-q');
    accuracySvm = sum(classesS == labs)/length(labs);
    
end

