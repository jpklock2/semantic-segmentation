function [labelImageTree] = runMasseli(modelSvm, modelTree, usefulImage)

    %% Feature Extraction

    % pegando retangulo do meio
%     x0 = max(utilCropSize(2), 1);
%     y0 = max(utilCropSize(1), 1);
%     x1 = min(utilCropSize(4), size(mask,2));
%     y1 = min(utilCropSize(3), size(mask,1));
% 
%     rec = [y0 x0 y1 x1];
%     plotMask(mask, classes, rec);
% 
%     usefulImage = rotPlotImage(y0:y1, x0:x1, :);
    usefulGray = rgb2gray(usefulImage);
    pad = 0;
    usefulGrayPad = padarray(usefulGray, [pad pad], 0, 'both');
    
    y1 = size(usefulGray, 1);
    x1 = size(usefulGray, 2);
    y0 = 1;
    x0 = 1;
    
    tY = y1-y0+2*pad;
    tX = x1-x0+2*pad;
    delta =170;
%     imSize = pad+pad+10;
    
%     figure;
%     imshow(usefulGray);
%     hold on;
    
    fts = [];
    labs = [];
    bestPoints = {[]};
    labelImageSvm = zeros(size(usefulGray));
    labelImageTree = zeros(size(usefulGray));
    
    totalY = length(pad+1:delta:tY-delta-pad);
    totalX = length(pad+1:delta:tX-delta-pad);
    
    cntY = 1;
    for ponY = pad+1:delta:tY-delta-pad
        fprintf('\nLinha %d/%d\n', cntY, totalY);
        cntX = 1;
        for ponX = pad+1:delta:tX-delta-pad
            fprintf('Coluna %d/%d\n', cntX, totalX);
            currImg = usefulGrayPad(ponY-pad:ponY-1+delta+pad, ponX-pad:ponX-1+delta+pad);
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
            bestPoints = [bestPoints; {bestP}];
            
            tempFts = extractFeatures(currImg, points);
%             currFts = double(tempFts.Features(bestIdx, :)); %./255;
            currFts = tempFts(bestIdx, :);
            fts = [fts; currFts];
            
%             labSvm = svmpredict(rand([1 1]), currFts, modelSvm, '-q');
            labTree = str2num(cell2mat(predict(modelTree, currFts)));
%             labTree = predict(modelTree, currFts);
%             labelImageSvm(ponY-pad:ponY-1+delta-pad, ponX-pad:ponX-1+delta-pad) = labSvm;
            labelImageTree(ponY-pad:ponY-1+delta-pad, ponX-pad:ponX-1+delta-pad) = labTree;
            
%             plot(bestP);
            cntX = cntX+1;
        end
        cntY = cntY+1;
    end
    
    figure; imagesc(labelImageSvm);
    figure; imagesc(labelImageTree);
    dbg = 1;
%     J = insertMarker(I,corners,'circle');
%     imshow(J)
     
end

