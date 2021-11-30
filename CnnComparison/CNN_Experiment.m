 %% Starting Code
initCode;
addpath('../SemanticSegmentation');
dataset = 'BR';
loadImageNames;
imageNamesTemp = imageNames(2:end);
loadColors;
classes = myColorsClasses';
cmap = myColors;

myColorsClasses = ["Grass", "Forest", "Soil", "River"];
myColorsNames = myColorsNames([2 4 8 1]);
myColors = myColors([2 4 8 1], :);
classes = myColorsClasses';
cmap = myColors;

trainImage = imread('CnnComparison/Images/Train/map.tif');
trainImage = imresize(trainImage, 0.25);
origTrainImage = trainImage;

labImage = rgb2lab(trainImage);
% [counts,binLocations] = imhist(rgbImage);
% figure; bar(binLocations,counts);
L = labImage(:,:,1)/100;
L = adapthisteq(L);
labImage(:,:,1) = L*100;
trainImage = lab2rgb(labImage);

% Compute histogram of the reference image
hgram = zeros(3, 256);
for currentChan = 1:3
    hgram(currentChan, :) = imhist(trainImage(:,:,currentChan), 256);
end

% [L, idx, centroidsFinal, classes, parameters, superPixels, pixelsOwn, pixelsAdj, rgbImage, originalRgbImage] = semanticSegmentation(imagePath, m, myFilter2, dataset, classes, centroidsFinal, parameters, scaleSize);

segmentMap = 0;
if segmentMap
[maskGeo, maskIdxGeo, centroidsGeo, classesGeo, parameters, geoAdjacencies, ftGeoOwn, ftGeoAdj, ~, currImagePlot] = semanticSegmentation('CnnComparison/Images/Train/map.tif', 1, 3, dataset, 0, 2);
LmaskGeo = zeros(size(maskGeo));
for mx = 1:length(classesGeo)
%     LmaskGeo(maskGeo == mx) = classesGeo(mx);
    LmaskGeo(maskIdxGeo{mx}) = classesGeo(mx);
end
% outputSegmentation = evalFunction(classesGeo, length(unique(classesGeo)), maskIdxGeo, currImagePlot, length(classesGeo));
save(['mapData.mat'],'maskGeo','maskIdxGeo','centroidsGeo','classesGeo','parameters','LmaskGeo');
else
    load('mapData.mat');
%     load('mapData_CANFIS.mat');
%     load('mapData_CANFIS_2.mat');
end

% [mask, maskIdx, centroids, classes, ~, adjacencies, ftOwn, ftAdj, currImage, currImagePlot] = semanticSegmentation(outf, i+1, 3, classesGeo, centroidsGeo, parameters, [360 480]);
% outputSegmentation = evalFunction(classes, length(unique(classes)), maskIdx, currImagePlot, length(classes));
% Lmask = zeros(size(mask));
% for mx = 1:length(classes)
%     Lmask(mask == mx) = classes(mx);
% end

%% Getting training data
generateDataset = 0;
if generateDataset
    
% trainImage = imread('Images/Train/Mosaicoo.tif');
% trainImage = imresize(trainImage, 0.25);
load('CnnComparison/Images/Train/Classes/myClasses.mat');

% combining classes
myClasses{1} = [myClasses{1}; myClasses{2}; myClasses{3}];
myClasses{2} = [myClasses{4}; myClasses{5}; myClasses{12}];
myClasses{3} = [myClasses{7}; myClasses{8}; myClasses{11}; myClasses{10}; myClasses{9}];
myClasses{4} = [myClasses{6}; myClasses{13}; myClasses{14}];
myClasses = myClasses(1:4);
% myColorsClasses = ["Grass", "Forest", "Soil", "River"];
% myColorsNames = myColorsNames([2 4 8 1]);
% myColors = myColors([2 4 8 1], :);
% classes = myColorsClasses;
% cmap = myColors;

% getting masks
idx = label2idx(L2);
Lmask = zeros(size(L2));
LmaskRed = zeros(size(L2));
LmaskGreen = zeros(size(L2));
LmaskBlue = zeros(size(L2));
LmaskColor = zeros(size(L2, 1), size(L2, 2), 3);
LmaskCategorical = strings(size(L2));
coloredImageRed = origTrainImage(:, :, 1);
coloredImageGreen = origTrainImage(:, :, 2);
coloredImageBlue = origTrainImage(:, :, 3);
coloredImage = origTrainImage;
sizesClasses = zeros(1, length(myClasses));
% numRows = size(L2, 1);
% numCols = size(L2, 2);
for mx = 1:length(myClasses)
    Lmask(ismember(L2, myClasses{mx})) = mx;
    LmaskRed(ismember(L2, myClasses{mx})) = myColors(mx, 1);
    LmaskGreen(ismember(L2, myClasses{mx})) = myColors(mx, 2);
    LmaskBlue(ismember(L2, myClasses{mx})) = myColors(mx, 3);
    LmaskCategorical(ismember(L2, myClasses{mx})) = myColorsClasses(mx);
    sizesClasses(mx) = length(myClasses{mx});
    coloredImageRed(ismember(L2, myClasses{mx})) = mean(coloredImageRed(ismember(L2, myClasses{mx})));
    coloredImageGreen(ismember(L2, myClasses{mx})) = mean(coloredImageGreen(ismember(L2, myClasses{mx})));
    coloredImageBlue(ismember(L2, myClasses{mx})) = mean(coloredImageBlue(ismember(L2, myClasses{mx})));
end
LmaskColor(:, :, 1) = LmaskRed;
LmaskColor(:, :, 2) = LmaskGreen;
LmaskColor(:, :, 3) = LmaskBlue;
coloredImage(:, :, 1) = coloredImageRed;
coloredImage(:, :, 2) = coloredImageGreen;
coloredImage(:, :, 3) = coloredImageBlue;
LmaskCategorical = categorical(LmaskCategorical, myColorsClasses);
% LmaskColor = [LmaskRed LmaskGren LmaskBlue];
% figure; imagesc(Lmask);
% figure; imshow(LmaskColor);

% B = labeloverlay(trainImage,LmaskCategorical,'Colormap',cmap,'Transparency',0.6);
% figure; imshow(B);
% pixelLabelColorbar(myColors,myColorsClasses);

% fig = figure; montage({origTrainImage, coloredImage, LmaskColor}, 'ThumbnailSize', [], 'BorderSize', [10 10], 'Size', [1 3], 'BackgroundColor', 'black');
% set(fig,'Units','Inches');
% pos = get(fig,'Position');
% set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(fig,'D:\2021_1\mestrado\artigo\SIBIGRAPI\images\dataset_BR.pdf','-dpdf','-r0');
% title(['Image ' num2str(i)]);

% Extracting superpixels centroids and adjacency matrix
numRows = size(L2,1);
numCols = size(L2,2);
superPixels = cell(1, N2);
for labelVal = 1:N2
    [y, x] = find(L2 == labelVal);
    superPixels{labelVal}.centX = round(mean(x));
    superPixels{labelVal}.centY = round(mean(y));
    superPixels{labelVal}.sizeX = max(x)-min(x);
    superPixels{labelVal}.sizeY = max(y)-min(y);
end

fprintf('\n');
imgCnt = 1;

% for i = 1:length(myClasses)
%     currClass = myClasses{i};
        currSamples = randsample(1:N2, 1000, true);
%     currSamples = randsample(currClass, 100, true);
    for j = 1:length(currSamples)
        fprintf('Image nº %d\n', imgCnt);
        samp = currSamples(j);
        
        
        % getting rectangle
        x0 = superPixels{samp}.centX-240+1;
        x1 = superPixels{samp}.centX+240;
        y0 = superPixels{samp}.centY-180+1;
        y1 = superPixels{samp}.centY+180;
        
        % assuring size
        if y0 < 1; y0 = 1; y1 = 360; elseif y1 > numRows; y0 = numRows-360+1; y1 = numRows; end
        if x0 < 1; x0 = 1; x1 = 480; elseif x1 > numCols; x0 = numCols-480+1; x1 = numCols; end
        
%         disp([y0 y1 x0 x1]);
%         disp(size(Lmask(y0:y1, x0:x1)));
        
        % getting image
        currImg = trainImage(y0:y1, x0:x1, :);
        currLabel = Lmask(y0:y1, x0:x1);
        currLabelColor = LmaskColor(y0:y1, x0:x1, :);
%         figure; imshow(currImg);
%         figure; imagesc(currLabel);
%         figure; imshow(currLabelColor);
%         currPath = '';
%         imwrite(uint8(currImg), ['CNNCodes/myDataset/images2/' num2str(imgCnt, '%05.f') '.png']);
        imwrite(currImg, ['CnnComparison/Images/CnnDataset/ImagesTrain/' num2str(imgCnt, '%05.f') '.png']);
        imwrite(currLabelColor, ['CnnComparison/Images/CnnDataset/LabelsTrain/' num2str(imgCnt, '%05.f') '.png']);
        imgCnt = imgCnt + 1;
        dbg = 1;
        
    end
% end

end

%% Reading images

imds = imageDatastore('CnnComparison/Images/CnnDataset/ImagesTrain');
% I = readimage(imds,25);
% I = histeq(I);
% imshow(I);

pxds = pixelLabelDatastore('CnnComparison/Images/CnnDataset/LabelsTrain', myColorsClasses, myColors*255);

% C = readimage(pxds,25);
% B = labeloverlay(I,C,'ColorMap',myColors);
% imshow(B);
% pixelLabelColorbar(myColors,myColorsClasses);

tbl = countEachLabel(pxds);
% frequency = tbl.PixelCount/sum(tbl.PixelCount);

% bar(1:numel(myColorsClasses),frequency);
% xticks(1:numel(myColorsClasses)) ;
% xticklabels(tbl.Name);
% xtickangle(45);
% ylabel('Frequency');

%% Convert Dataset

convertDataset = 0;
if convertDataset
for m = 1:10
    
currentImage = strtrim(imageNamesTemp{m});   
rgbImage = imread(currentImage);
    
labImage = rgb2lab(rgbImage);
% [counts,binLocations] = imhist(rgbImage);
% figure; bar(binLocations,counts);
L = labImage(:,:,1)/100;
L = adapthisteq(L);
labImage(:,:,1) = L*100;
rgbImage = lab2rgb(labImage);

for k = 1:size(rgbImage, 3) % Process one color channel at a time
    hgramToUse = k;
    rgbImage(:, :, k) = histeq(rgbImage(:,:,k), hgram(hgramToUse,:));
end

rgbImage = imresize(rgbImage,[360, 480]);

load(['CnnComparison/Images/Test/Original_Dev/Classes/im' num2str(m+1) '.mat']);

% combining classes
myClasses{1} = [myClasses{1}; myClasses{2}; myClasses{3}];
myClasses{2} = [myClasses{4}; myClasses{5}; myClasses{12}];
myClasses{3} = [myClasses{7}; myClasses{8}; myClasses{11}; myClasses{10}; myClasses{9}];
myClasses{4} = [myClasses{6}; myClasses{13}; myClasses{14}];
myClasses = myClasses(1:4);
% myColorsClasses = ["Grass", "Forest", "Soil", "River"];
% myColorsNames = myColorsNames([2 4 8 1]);
% myColors = myColors([2 4 8 1], :);
% classes = myColorsClasses;
% cmap = myColors;

% getting masks
idx = label2idx(L2);
Lmask = zeros(size(L2));
LmaskRed = zeros(size(L2));
LmaskGreen = zeros(size(L2));
LmaskBlue = zeros(size(L2));
LmaskColor = zeros(size(L2, 1), size(L2, 2), 3);
% LmaskCategorical = strings(size(L2));
sizesClasses = zeros(1, length(myClasses));
for mx = 1:length(myClasses)
    Lmask(ismember(L2, myClasses{mx})) = mx;
    LmaskRed(ismember(L2, myClasses{mx})) = myColors(mx, 1);
    LmaskGreen(ismember(L2, myClasses{mx})) = myColors(mx, 2);
    LmaskBlue(ismember(L2, myClasses{mx})) = myColors(mx, 3);
%     LmaskCategorical(ismember(L2, myClasses{mx})) = myColorsClasses(mx);
    sizesClasses(mx) = length(myClasses{mx});
end
LmaskColor(:, :, 1) = LmaskRed;
LmaskColor(:, :, 2) = LmaskGreen;
LmaskColor(:, :, 3) = LmaskBlue;
% LmaskCategorical = categorical(LmaskCategorical, myColorsClasses);
LmaskColorTemp = imresize(LmaskColor, [360, 480], 'nearest');

% LmaskCategorical = strings([360, 480]);
% for mx = 1:length(myColors)
%     members = ismember(LmaskColorTemp, myColors(mx, :));
%     LmaskCategorical(members(:, :, 1)) = myColorsClasses(mx);
% end
% LmaskCategorical = categorical(LmaskCategorical, myColorsClasses);
% expectedResult = LmaskCategorical;

imwrite(rgbImage, ['CnnComparison/Images/CnnDataset/ImagesTest/' num2str(m, '%05.f') '.png']);
imwrite(LmaskColorTemp, ['CnnComparison/Images/CnnDataset/LabelsTest/' num2str(m, '%05.f') '.png']);
dbg = 1;

end
end

%% Reading test dataset
% currentImage = strtrim(imageNamesTemp{m});   
% rgbImage = imread(currentImage);
% originalRgbImage = rgbImage;
% 
% load(['Images/Test/Original_Dev/Classes/im' num2str(m+1) '.mat']);
% 
% % combining classes
% myClasses{1} = [myClasses{1}; myClasses{2}; myClasses{3}];
% myClasses{2} = [myClasses{4}; myClasses{5}; myClasses{12}];
% myClasses{3} = [myClasses{7}; myClasses{8}; myClasses{11}; myClasses{10}; myClasses{9}];
% myClasses{4} = [myClasses{6}; myClasses{13}; myClasses{14}];
% myClasses = myClasses(1:4);
% myColorsClasses = ["Grass", "Forest", "Soil", "River"];
% myColorsNames = myColorsNames([2 4 8 1]);
% myColors = myColors([2 4 8 1], :);
% classes = myColorsClasses;
% cmap = myColors;

% getting masks
% idx = label2idx(L2);
% Lmask = zeros(size(L2));
% LmaskRed = zeros(size(L2));
% LmaskGren = zeros(size(L2));
% LmaskBlue = zeros(size(L2));
% LmaskColor = zeros(size(L2, 1), size(L2, 2), 3);
% % LmaskCategorical = strings(size(L2));
% sizesClasses = zeros(1, length(myClasses));
% for mx = 1:length(myClasses)
%     Lmask(ismember(L2, myClasses{mx})) = mx;
%     LmaskRed(ismember(L2, myClasses{mx})) = myColors(mx, 1);
%     LmaskGren(ismember(L2, myClasses{mx})) = myColors(mx, 2);
%     LmaskBlue(ismember(L2, myClasses{mx})) = myColors(mx, 3);
% %     LmaskCategorical(ismember(L2, myClasses{mx})) = myColorsClasses(mx);
%     sizesClasses(mx) = length(myClasses{mx});
% end
% LmaskColor(:, :, 1) = LmaskRed;
% LmaskColor(:, :, 2) = LmaskGren;
% LmaskColor(:, :, 3) = LmaskBlue;
% LmaskCategorical = categorical(LmaskCategorical, myColorsClasses);

% LmaskColorTemp = imresize(LmaskColor, [360, 480], 'nearest');
% LmaskCategorical = strings([360, 480]);
% for mx = 1:length(myColors)
%     members = ismember(LmaskColorTemp, myColors(mx, :));
%     LmaskCategorical(members(:, :, 1)) = myColorsClasses(mx);
% end
% LmaskCategorical = categorical(LmaskCategorical, myColorsClasses);
% expectedResult = LmaskCategorical;

% I = histeq(originalRgbImage);

imdsTest = imageDatastore('CnnComparison/Images/CnnDataset/ImagesTest');

pxdsTest = pixelLabelDatastore('CnnComparison/Images/CnnDataset/LabelsTest', myColorsClasses, myColors*255);

% labImage = rgb2lab(originalRgbImage);
% % [counts,binLocations] = imhist(rgbImage);
% % figure; bar(binLocations,counts);
% L = labImage(:,:,1)/100;
% L = adapthisteq(L);
% labImage(:,:,1) = L*100;
% I2 = lab2rgb(labImage);

% figure; montage({originalRgbImage, I, I2});

% for k = 1:size(rgbImage, 3) % Process one color channel at a time
%     hgramToUse = k;
%     rgbImage(:, :, k) = histeq(rgbImage(:,:,k), hgram(hgramToUse,:));
% end
% figure; montage({originalRgbImage, rgbImage});

%% Generate FCM Dataset (no Classifier)

generateFCM_ori = 0;
if generateFCM_ori
    
listing = dir('CnnComparison/Images/Test/Original_Dev');
imageNames = [{}];
for i=1:length(listing)
    if contains(lower(listing(i).name), '.jpg')
        imageNames = [imageNames; listing(i).name];
    end
end

currImagePlots = [{}];
outputImages = [{}];
LmaskColors = [{}];
for m = 1:length(imageNames)

    [~, maskIdx, ~, classes, params, ~, ~, ~, ~, currImagePlot] = semanticSegmentation(imageNames{m}, 1, 3);
    
    outputImage = zeros(size(currImagePlot), 'like', currImagePlot);
    numRows = size(outputImage, 1);
    numCols = size(outputImage, 2);
    for labelVal = 1:length(maskIdx)
        redIdx = maskIdx{labelVal};
        greenIdx = maskIdx{labelVal}+numRows*numCols;
        blueIdx = maskIdx{labelVal}+2*numRows*numCols;
        outputImage(redIdx) = params.meanRed(classes(labelVal));
        outputImage(greenIdx) = params.meanGreen(classes(labelVal));
        outputImage(blueIdx) = params.meanBlue(classes(labelVal));
    end
    
    LmaskColor = imread(pxdsTest.Files{m});
    LmaskColor = imresize(LmaskColor, [size(currImagePlot, 1), size(currImagePlot, 2)], 'Bilinear');
    
    currImagePlots = [currImagePlots; currImagePlot];
    outputImages = [outputImages; outputImage];
    LmaskColors = [LmaskColors; LmaskColor];
    imwrite(LmaskColor, ['CnnComparison/Images/Results/Ours/' num2str(m, '%05.f') '.png']);
%     fig = figure;
%     montage({currImagePlot, outputImage, LmaskColor}, 'ThumbnailSize', [], 'BorderSize', [5 5], 'Size', [1 3], 'BackgroundColor', 'black');
%     set(fig,'Units','Inches');
%     pos = get(fig,'Position');
%     set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(fig,['D:\2021_1\mestrado\artigo\SIBIGRAPI\images\ours_' num2str(m) '.pdf'],'-dpdf','-r0');
%     title(['Image ' num2str(i)]);

end
end

% fig = figure; montage({currImagePlots{1}, outputImages{1}, LmaskColors{1}, currImagePlots{2}, outputImages{2}, LmaskColors{2},...
%         currImagePlots{4}, outputImages{4}, LmaskColors{4}},...
%     'ThumbnailSize', [], 'BorderSize', [5 5], 'Size', [3 3], 'BackgroundColor', 'black');
%     currImagePlots{4}, outputImages{4}, LmaskColors{4}, currImagePlots{7}, outputImages{7}, LmaskColors{7}},...
%     'ThumbnailSize', [], 'BorderSize', [5 5], 'Size', [4 3], 'BackgroundColor', 'black');
%     set(fig,'Units','Inches');
%     pos = get(fig,'Position');
%     set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(fig,'D:\2021_1\mestrado\artigo\SIBIGRAPI\images\sem_segs.pdf','-dpdf','-r0');

% fprintf("\nOurs Reslts...\n");
% pxdsOurs = pixelLabelDatastore('CnnComparison/Images/Results/Ours', myColorsClasses, myColors*255);
% metrics = evaluateSemanticSegmentation(pxdsOurs,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics

% fprintf("\nCANFIS 2 Centroids No Dendo Reslts...\n");
% pxdsFCM = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_CANFIS_2_centroids_nodendo', myColorsClasses, myColors*255);
% metrics = evaluateSemanticSegmentation(pxdsFCM,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics

% pxdsOurs = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_SVM_superpixels', myColorsClasses, myColors*255);

dbg = 1;

%% Generate FCM Dataset

generateFCM = 0;
if generateFCM
for m = 1:10

    [mask, maskIdx, centroids, classes, ~, adjacencies, ftOwn, ftAdj, currImage, currImagePlot] = semanticSegmentation(imdsTest.Files{m}, m+1, 3, '', 1, 2, classesGeo, centroidsGeo, parameters, [360 480]);
    
    %     Lmask = zeros(size(mask));
    LmaskRed = zeros(size(mask));
    LmaskGreen = zeros(size(mask));
    LmaskBlue = zeros(size(mask));
    LmaskColor = zeros(size(mask, 1), size(mask, 2), 3);
%     LmaskCategorical = strings(size(mask));
    %     myClasses = cell(1, 4);
    for mx = 1:length(classes)
    %         Lmask(mask == mx) = classes(mx);
    %         myClasses{classes(mx)} = [myClasses{classes(mx)}; maskIdx{mx}];
%         LmaskRed(mask == mx) = myColors(classes(mx), 1);
%         LmaskGreen(mask == mx) = myColors(classes(mx), 2);
%         LmaskBlue(mask == mx) = myColors(classes(mx), 3);
        
        LmaskRed(maskIdx{mx}) = myColors(classes(mx), 1);
        LmaskGreen(maskIdx{mx}) = myColors(classes(mx), 2);
        LmaskBlue(maskIdx{mx}) = myColors(classes(mx), 3);
%         LmaskCategorical(mask == mx) = myColorsClasses(classes(mx));
    end
    LmaskColor(:, :, 1) = LmaskRed;
    LmaskColor(:, :, 2) = LmaskGreen;
    LmaskColor(:, :, 3) = LmaskBlue;
%     LmaskCategorical = categorical(LmaskCategorical, myColorsClasses);
    
    imwrite(LmaskColor, ['CnnComparison/Images/Results/Ours/' num2str(m, '%05.f') '.png']);
    
end
end

fprintf("\nOurs Reslts...\n");
pxdsOurs = pixelLabelDatastore('CnnComparison/Images/Results/Ours', myColorsClasses, myColors*255);
metrics = evaluateSemanticSegmentation(pxdsOurs,pxdsTest,'Verbose',false);
metrics.DataSetMetrics
metrics.ClassMetrics

% fprintf("\nSVM Superpixels Reslts...\n");
% pxdsFCM = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_SVM_superpixels', myColorsClasses, myColors*255);
% metrics = evaluateSemanticSegmentation(pxdsFCM,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics
% 
% fprintf("\nSVM Superpixels Reslts No Increased K...\n");
% pxdsFCM = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_SVM_superpixels_noincreased', myColorsClasses, myColors*255);
% metrics = evaluateSemanticSegmentation(pxdsFCM,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics
% 
% fprintf("\nSVM Centroids No Dendo Reslts...\n");
% pxdsFCM = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_SVM_centroids_nodendo', myColorsClasses, myColors*255);
% metrics = evaluateSemanticSegmentation(pxdsFCM,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics
% 
% fprintf("\nSVM Centroids With Dendo Reslts...\n");
% pxdsFCM = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_SVM_centroids_withdendo', myColorsClasses, myColors*255);
% metrics = evaluateSemanticSegmentation(pxdsFCM,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics
% 
% fprintf("\nCANFIS Centroids No Dendo Reslts...\n");
% pxdsFCM = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_CANFIS_centroids_nodendo', myColorsClasses, myColors*255);
% metrics = evaluateSemanticSegmentation(pxdsFCM,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics
% % 
% fprintf("\nCANFIS Superpixels Reslts...\n");
% pxdsFCM = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_CANFIS_superpixels', myColorsClasses, myColors*255);
% metrics = evaluateSemanticSegmentation(pxdsFCM,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics
% 
% fprintf("\nCANFIS 2 Superpixels Reslts...\n");
% pxdsFCM = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_CANFIS_2_superpixels', myColorsClasses, myColors*255);
% metrics = evaluateSemanticSegmentation(pxdsFCM,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics
% 
% fprintf("\nCANFIS 2 Centroids No Dendo Reslts...\n");
% pxdsFCM = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_CANFIS_2_centroids_nodendo', myColorsClasses, myColors*255);
% metrics = evaluateSemanticSegmentation(pxdsFCM,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics

% pxdsOurs = pixelLabelDatastore('CNNCodes/myDataset/labelsFCM_SVM_superpixels', myColorsClasses, myColors*255);

dbg = 1;

%% FCN

trainFCN = 0;
if trainFCN
imageSize = [360 480];
numClasses = numel(myColorsClasses);
lgraph = fcnLayers(imageSize,numClasses);

imageFreq = tbl.PixelCount ./ tbl.ImagePixelCount;
classWeights = median(imageFreq) ./ imageFreq;

% pxLayer = pixelClassificationLayer('Name','labels','Classes',tbl.Name,'ClassWeights',classWeights);
pxLayer = pixelClassificationLayer('Name','labels','Classes',tbl.Name);

lgraph = removeLayers(lgraph,'pixelLabels');
lgraph = addLayers(lgraph, pxLayer);
lgraph = connectLayers(lgraph,'softmax','labels');

options = trainingOptions('adam', ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',1/3, ...
    'LearnRateDropPeriod',5, ...
    'MaxEpochs',50, ...  
    'MiniBatchSize',3, ...
    'Shuffle','every-epoch', ...
    'CheckpointPath', 'D:/temp', ...
    'VerboseFrequency',1, ...
    'Plots','training-progress');

% augmenter = imageDataAugmenter('RandXReflection',true,...
%     'RandXTranslation',[-10 10],'RandYTranslation',[-10 10],...
%     'RandRotation', [-45 45], 'RandScale', [0.5 4]);

pximds = pixelLabelImageDatastore(imds,pxds); %, ...
%     'DataAugmentation',augmenter);

doTraining = true;
if doTraining
%     load('D:\temp\net_checkpoint__2997__2021_07_07__03_56_51.mat','net');
%     newlgraph = layerGraph(net);
%     [net2, info2] = trainNetwork(pximds,newlgraph,options);
    [net2, info2] = trainNetwork(pximds,lgraph,options);
    save('CnnComparison/TrainedNetworks/fcn.mat','net2','info2');
end
end

%% DeepLab

trainDeepLab = 0;
if trainDeepLab
% Specify the network image size. This is typically the same as the traing image sizes.
imageSize = [360 480 3];

% Specify the number of classes.
numClasses = numel(classes);

% Create DeepLab v3+.
lgraph = deeplabv3plusLayers(imageSize, numClasses, "inceptionresnetv2");
pxLayer = pixelClassificationLayer('Name','labels','Classes',tbl.Name);
lgraph = replaceLayer(lgraph,"classification",pxLayer);

% Define validation data.
% dsVal = combine(imdsVal,pxdsVal);

% Define training options. 
options = trainingOptions('adam', ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',1/3, ...
    'LearnRateDropPeriod',5, ...
    'MaxEpochs',50, ...  
    'MiniBatchSize',3, ...
    'Shuffle','every-epoch', ...
    'CheckpointPath', 'D:/temp', ...
    'VerboseFrequency',1, ...
    'Plots','training-progress');

dsTrain = combine(imds, pxds);

doTraining = true;
if doTraining    
    [net, info] = trainNetwork(dsTrain,lgraph,options);
    save('CnnComparison/TrainedNetworks/deeplab.mat','net','info');
else
    data = load(pretrainedNetwork); 
    net = data.net;
end

end

%% SegNet

trainSegNet = 0;
if trainSegNet
% Specify the network image size. This is typically the same as the traing image sizes.
imageSize = [360 480 3];

% Specify the number of classes.
numClasses = numel(myColorsClasses);

% Create SegNet
lgraph = segnetLayers(imageSize,numClasses,'vgg19');
pxLayer = pixelClassificationLayer('Name','labels','Classes',tbl.Name);
lgraph = removeLayers(lgraph,'pixelLabels');
lgraph = addLayers(lgraph, pxLayer);
lgraph = connectLayers(lgraph,'softmax','labels');

% Define validation data.
% dsVal = combine(imdsVal,pxdsVal);

% Define training options. 
options = trainingOptions('adam', ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',1/3, ...
    'LearnRateDropPeriod',5, ...
    'MaxEpochs',50, ...  
    'MiniBatchSize',3, ...
    'Shuffle','every-epoch', ...
    'CheckpointPath', 'D:/temp', ...
    'VerboseFrequency',1, ...
    'Plots','training-progress');

dsTrain = combine(imds, pxds);

doTraining = true;
if doTraining    
    [net, info] = trainNetwork(dsTrain,lgraph,options);
    save('CnnComparison/TrainedNetworks/segnet.mat','net','info');
else
    data = load(pretrainedNetwork); 
    net = data.net;
end

end

%% SegNet
runSegNet = 1;
if runSegNet
% rgbImage = imresize(originalRgbImage,[360, 480]);
% rgbImage = im2double(rgbImage);
% imshow(rgbImage);
% predict_scores = segnet_predict_mex(rgbImage);
% [~,argmax] = max(predict_scores,[],3);
% classes = [
%     "Sky"
%     "Building"
%     "Pole"
%     "Road"
%     "Pavement"
%     "Tree"
%     "SignSymbol"
%     "Fence"
%     "Car"
%     "Pedestrian"
%     "Bicyclist"
%     ];
% cmap = camvidColorMap();
% SegmentedImage = labeloverlay(rgbImage,argmax,'ColorMap',cmap);
% figure;
% imshow(SegmentedImage);
% pixelLabelColorbar(cmap,classes);
% title('SegNet');

fprintf("\nSegNet Reslts...\n");
data = load('CnnComparison/TrainedNetworks/segnet.mat'); 
segnet = data.net;
segResults = semanticseg(imdsTest,segnet, ...
    'MiniBatchSize',3, ...
    'WriteLocation','D:/temp', ...
    'Verbose',false);
metrics = evaluateSemanticSegmentation(segResults,pxdsTest,'Verbose',false);
metrics.DataSetMetrics
metrics.ClassMetrics

end

dbg = 1;

%% ResNet
runResNet = 1;
if runResNet
% pretrainedURL = 'https://www.mathworks.com/supportfiles/vision/data/deeplabv3plusResnet18CamVid.mat';
% pretrainedFolder = fullfile('CNNCodes');
% pretrainedNetwork = fullfile(pretrainedFolder,'deeplab_satellite.mat'); 
% if ~exist(pretrainedNetwork,'file')
%     mkdir(pretrainedFolder);
%     disp('Downloading pretrained network (58 MB)...');
%     websave(pretrainedNetwork,pretrainedURL);
% end
% Specify the network image size. This is typically the same as the traing image sizes.
% imageSize = [720 960 3];
% Specify the number of classes.
% numClasses = numel(classes);
% Create DeepLab v3+.
% lgraph = deeplabv3plusLayers(imageSize, numClasses, "resnet18");
% imageFreq = tbl.PixelCount ./ tbl.ImagePixelCount;
% classWeights = median(imageFreq) ./ imageFreq;
% pxLayer = pixelClassificationLayer('Name','labels','Classes',classes);%,'ClassWeights',classWeights);
% lgraph = replaceLayer(lgraph,"classification",pxLayer);
% Define validation data.
% dsVal = combine(imdsVal,pxdsVal);
% Define training options. 
% options = trainingOptions('sgdm', ...
%     'LearnRateSchedule','piecewise',...
%     'LearnRateDropPeriod',10,...
%     'LearnRateDropFactor',0.3,...
%     'Momentum',0.9, ...
%     'InitialLearnRate',1e-3, ...
%     'L2Regularization',0.005, ...
%     %'ValidationData',dsVal,...;
%     'MaxEpochs',30, ...  
%     'MiniBatchSize',8, ...
%     'Shuffle','every-epoch', ...
%     'CheckpointPath', 'D:/temp', ...
%     'VerboseFrequency',2,...
%     'Plots','training-progress',...
%     'ValidationPatience', 4);
% dsTrain = combine(imdsTrain, pxdsTrain);
% xTrans = [-10 10];
% yTrans = [-10 10];
% dsTrain = transform(dsTrain, @(data)augmentImageAndLabel(data,xTrans,yTrans));
% doTraining = false;
% if doTraining    
%     [net, info] = trainNetwork(dsTrain,lgraph,options);
% else
%     data = load(pretrainedNetwork); 
%     net = data.net;
% end
% rgbImage = imresize(rgbImage,[360 480]);
% rgbImage = im2double(rgbImage);
% imshow(rgbImage);
% I = readimage(imdsTest,35);
% C = semanticseg(rgbImage, net);
% B = labeloverlay(rgbImage,C,'Colormap',cmap,'Transparency',0.4);
% figure;
% imshow(B)
% pixelLabelColorbar(cmap, classes);
% title('ResNet');
% expectedResult = readimage(pxdsTest,35);
% actual = uint8(C);
% expected = uint8(expectedResult);
% imshowpair(actual, expected)
% iou = jaccard(C,expectedResult);
% table(classes,iou)
% pxdsResults = semanticseg(imdsTest,net, ...
%     'MiniBatchSize',4, ...
%     'WriteLocation',tempdir, ...
%     'Verbose',false);
% metrics = evaluateSemanticSegmentation(pxdsResults,pxdsTest,'Verbose',false);
% metrics.DataSetMetrics
% metrics.ClassMetrics

fprintf("\nDeepLab Reslts...\n");
data = load('CnnComparison/TrainedNetworks/deeplab.mat'); 
deepnet = data.net;
deepResults = semanticseg(imdsTest,deepnet, ...
    'MiniBatchSize',3, ...
    'WriteLocation','D:/temp', ...
    'Verbose',false);
metrics = evaluateSemanticSegmentation(deepResults,pxdsTest,'Verbose',false);
metrics.DataSetMetrics
metrics.ClassMetrics

end

dbg = 1;

%% FCN
runFCN = 1;
if runFCN
% doTraining = false;
% if ~doTraining
%     pretrainedURL = 'https://www.mathworks.com/supportfiles/gpucoder/cnn_models/fcn/FCN8sCamVid.mat';
%     disp('Downloading pretrained FCN (448 MB)...');
%     websave('FCN8sCamVid.mat',pretrainedURL);
% end
% rgbImage = imresize(originalRgbImage,[360, 480]);

% labImage = rgb2lab(rgbImage);
% [counts,binLocations] = imhist(rgbImage);
% figure; bar(binLocations,counts);
% L = labImage(:,:,1)/100;
% L = adapthisteq(L);
% labImage(:,:,1) = L*100;
% rgbImage = lab2rgb(labImage);
% 
% for k = 1:size(rgbImage, 3) % Process one color channel at a time
%     hgramToUse = k;
%     rgbImage(:, :, k) = histeq(rgbImage(:,:,k), hgram(hgramToUse,:));
% end

fprintf("\nFCN Reslts...\n");

data = load('CnnComparison/TrainedNetworks/fcn.mat'); 
fcnnet = data.net2;
% predict_scores = fcn_predict_mex(rgbImage);
% I = readimage(imdsTest,1);
% C = semanticseg(I, net);
% B = labeloverlay(I,C,'Colormap',cmap,'Transparency',0.4);
% imshow(B)
% pixelLabelColorbar(cmap, classes);
% expectedResult = readimage(pxdsTest,1);
% actual = uint8(C);
% expected = uint8(expectedResult);
% imshowpair(actual, expected)
% iou = jaccard(C,expectedResult);
% table(classes,iou)

fcnResults = semanticseg(imdsTest,fcnnet, ...
    'MiniBatchSize',3, ...
    'WriteLocation','D:/temp', ...
    'Verbose',false);

metrics = evaluateSemanticSegmentation(fcnResults,pxdsTest,'Verbose',false);
metrics.DataSetMetrics
metrics.ClassMetrics

% [~,argmax] = max(predict_scores,[],3);
% classes = [
%     "Sky"
%     "Building"
%     "Pole"
%     "Road"
%     "Pavement"
%     "Tree"
%     "SignSymbol"
%     "Fence"
%     "Car"
%     "Pedestrian"
%     "Bicyclist"
%     ];

% cmap = camvidColorMap();
% SegmentedImage = labeloverlay(rgbImage,argmax,'ColorMap',cmap);
% figure;
% imshow(SegmentedImage);
% pixelLabelColorbar(cmap,classes);
% title('FCN');
end

%% Visually checking results

segTimes = [];
deepTimes = [];
fcnTimes = [];
fcmTimes = [];
% for j = 1:10
j = 1;
for i = 1:10
    
%     fprintf(['Image ' num2str(i) ' experiment ' num2str(j) '\n']);
    I = readimage(imdsTest,i);
    
%     tic;
%     [mask, maskIdx, centroids, classes, ~, adjacencies, ftOwn, ftAdj, currImage, currImagePlot] = semanticSegmentation(imdsTest.Files{i}, i+1, 3, '', classesGeo, centroidsGeo, parameters, [360 480]);
%     fcmTimes = [fcmTimes; toc];
    
%     LmaskCategorical = strings(size(mask));
%     for mx = 1:length(classes)
%         LmaskCategorical(mask == mx) = myColorsClasses(classes(mx));
%     end
%     LmaskCategorical = categorical(LmaskCategorical, myColorsClasses);
%     C0 = LmaskCategorical;
    C0 = readimage(pxdsOurs,i);
    B0 = labeloverlay(I,C0,'ColorMap',myColors,'Transparency',0.4);
    
    C1 = readimage(pxdsTest,i);
    B1 = labeloverlay(I,C1,'ColorMap',myColors,'Transparency',0.4);
    
%     tic;
    C2 = semanticseg(I, segnet, 'ExecutionEnvironment', 'cpu');
%     segTimes = [segTimes; toc];
    B2 = labeloverlay(I,C2,'Colormap',cmap,'Transparency',0.4);
    
%     tic;
    C3 = semanticseg(I, deepnet, 'ExecutionEnvironment', 'cpu');
%     deepTimes = [deepTimes; toc];
    B3 = labeloverlay(I,C3,'Colormap',cmap,'Transparency',0.4);
    
%     tic;
%     C4 = semanticseg(I, fcnnet, 'ExecutionEnvironment', 'cpu');
%     fcnTimes = [fcnTimes; toc];
%     B4 = labeloverlay(I,C4,'Colormap',cmap,'Transparency',0.4);
    
    fig = figure; montage({B1, B2, B3, B0}, 'BorderSize', [5 5], 'ThumbnailSize', []);
%     subplot(221); imshow(B2);
%     subplot(222); imshow(B3);
%     subplot(223); imshow(B4);
%     subplot(224); imshow(B0);
    pixelLabelColorbar(cmap,classes);
%     set(fig,'Units','Inches');
%     pos = get(fig,'Position');
%     set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%     print(fig,['D:\2021_1\mestrado\artigo\SIBIGRAPI\images\deep_' num2str(i) '.pdf'],'-dpdf','-r0');
    
%     print(fig,['D:\2021_1\mestrado\artigo\deep.pdf'],'-dpdf','-r0');
%     title(['Image ' num2str(i)]);

    dbg = 1;

end
% end

dbg = 1;

% end