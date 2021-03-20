if printResults
fprintf('\nDefining Cluster Number for Validation Image...\n');
end
tic;
jotas = [];
for k = 2:30
    jotasInd = [];
    for k2 = 1:10
        [~, ~, Jc, ~] = MyFuzzyMeans_opt(pixels, k);
        jotasInd = [jotasInd; Jc(end)];
    end
    jotas = [jotas; median(jotasInd)];
end
Kc = knee_pt(jotas);
if Kc <= 0
    Kc = 1;
end
% Kc = 7;
% for Kc = 20:-1:10
if printResults
fprintf('Cluster Number = %d\n', Kc);
end

Nexp = 1000;
allUc = cell(Nexp, 1);
allJc = zeros(Nexp, 1);
allCentroids = cell(Nexp, 1);
allMetrics = zeros(Nexp, 3);
for k = 1:Nexp
    [idxClust, allUc{k}, JcTemp, allCentroids{k}] = MyFuzzyMeans_opt(pixels, Kc);
    allJc(k) = JcTemp(end);
    allMetrics(k, 1) = evalclusters(pixels, idxClust, 'CalinskiHarabasz').CriterionValues; % maior, melhor
    allMetrics(k, 2) = evalclusters(pixels, idxClust, 'DaviesBouldin').CriterionValues; % menor, melhor
    allMetrics(k, 3) = evalclusters(pixels, idxClust, 'silhouette').CriterionValues; % maior, melhor
end
idxCal = find(allMetrics(:, 1) == max(allMetrics(:, 1)), 1);
idxDav = find(allMetrics(:, 2) == min(allMetrics(:, 2)), 1);
idxSil = find(allMetrics(:, 3) == max(allMetrics(:, 3)), 1);
medianJc = abs(allJc - median(allJc));
meanJc = abs(allJc - mean(allJc));
idxJcMedian = find(medianJc == min(medianJc), 1);
idxJcMean = find(meanJc == min(meanJc), 1);
% idxJc = idxJcMedian;
% idxJc = idxJcMean;
% idxJc = find(allJc == min(allJc), 1);
idxJc = idxCal;
% idxJc = idxDav;
% idxJc = idxSil;

methods = {'centroid', 'median', 'ward', 'average', 'complete', 'single', 'weighted'};
metrics = {'euclidean', 'cityblock', 'chebychev', 'cosine', 'correlation'};
cDist = zeros(length(methods), length(metrics));
% cClust = zeros(length(methods), length(metrics));
for meth = 1:length(methods)
    for metr = 1:length(metrics)
        if meth < 4 && metr > 1
            continue;
        end
        D = pdist(allCentroids{idxJc}, metrics{metr});
        tree = linkage(allCentroids{idxJc}, methods{meth}, metrics{metr});
        cDist(meth, metr) = cophenet(tree, D);
%         cClust(meth, metr) = cophenet(tree, allCentroids{idxJc});
    end
end

[maxMeth, maxMetr] = find(cDist == max(max(cDist)), 1);
% [maxMeth, maxMetr] = find(cClust == max(max(cClust)), 1);
% maxMeth = 4; maxMetr = 6;
D = pdist(allCentroids{idxJc}, metrics{maxMetr});
tree = linkage(allCentroids{idxJc}, methods{maxMeth}, metrics{maxMetr});

% leafOrder = optimalleaforder(tree, D);
% fig = figure;
% [H, T, OUTPERM] = dendrogram(tree, 'Reorder', leafOrder);
% 
% imagem 3 é boa!
% xlabel('Cluster Number', 'fontsize', 15);
% ylabel('Height', 'fontsize', 15);
% title('Dendrogram of a Clustering Solution ', 'fontsize', 20);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',15);
% hold on; thLine = plot([min(leafOrder)-1 max(leafOrder+1)], [0.3*max(tree(:, 3)) 0.3*max(tree(:, 3))], 'r', 'Linewidth', 3);
% legend([thLine],'30% of' + string(newline) + 'Tree Height');

dbg = 1;

% sizess = [10 5^2 10^2 size(pixels, 1)];
% patchess = [10 25 50 75 100];
% metrics = zeros(Nexp, length(sizess)*length(patchess));
% outputImages = [{}];
% for k = 1:Nexp
%     [~, idxUc] = sort(Uc{k}, 2);
%     classesTemp = idxUc(:, end);
%     if plotsCompare
%         outputImages = [outputImages; evalFunction(classesTemp, Kc, idx, rgbImage, N)];
% %         [L, k, nodeRows, nodeCols] = superpixelSpacing(outputImages{k}, N);
%         figure; montage({rgbImage, outputImages{k}});
%     end
% %     metrics(k, 1) = std(outputImages{k}(:));
% %     metrics(k, 2) = entropy(outputImage(:));
% %     glcms = graycomatrix(im2gray(outputImage));
% %     stats = graycoprops(glcms);
% %     metrics(k, 3) = stats.Contrast;
% %     metrics(k, 4) = stats.Correlation;
% %     metrics(k, 5) = stats.Energy;
% %     metrics(k, 6) = stats.Homogeneity;
% %     metrics(k, 7) = sum(outputImage(:).^2);
%     cntM = 1;
%     for k1 = 1:length(sizess)
%         for k2 = 1:length(patchess)
% %             disp([8+cntM sizess(k1) patchess(k2)]);
%             metrics(k, cntM) = calculateHomogeinity(outputImages{k}, sizess(k1), patchess(k2));
%             cntM = cntM+1;
%         end
%     end
% end
% 
% idxJc = find(metrics(:,7) == min(metrics(:,7)), 1);

Uc = allUc{idxJc};
centroids = allCentroids{idxJc};
[~, idxUc] = sort(Uc, 2);
classesTemp = idxUc(:, end);
[classesTempJoin, KcJoin, UcJoin, centroidsJoin] = joinClasses(tree, Kc, classesTemp, centroids, Uc, pixels);

% [outputImage] = evalFunction(classesTempJoin, KcJoin, idx, originalRgbImage, N);
% fig = figure; montage({outputImageLabel,outputImage}, 'Size', [1 2]);

if printResults
fprintf('Execution time for defining clusters: %f s\n', toc);
end
if plotsCompare
%     thresh = 0.15;
%     [outputImage] = evalFunction(classesTemp, Kc, idx, rgbImage, N);
    [outputImage] = evalFunction(classesTempJoin, KcJoin, idx, originalRgbImage, N);
    
%     load(['./Images/Test/Original_Dev/Classes/im',num2str(m),'.mat']);
%     idxLabel = label2idx(L2);
%     classesLabel = zeros(N2, 1);
%     for i = 1:length(myClasses)
%         currC = myClasses{i};
%         classesLabel(currC) = i;
%     end
%     [outputImageLabel] = evalFunction(classesLabel, 14, idxLabel, originalRgbImage, N2);
    
%     h9 = figure(5);
%     BW3 = boundarymask(L2);
%     h10 = imshow(imoverlay(rgbImage,BW3,'cyan'),'InitialMagnification',100);
%     rgbImageCrop2 = imoverlay(rgbImage,BW3,'cyan');
%     hold on;
%     h11 = imshow(tempImage);
%     h11.AlphaData = 0.75;
%     h12 = figure(6);
%     h13 = imshow(imoverlay(rgbImage,BW3,'cyan'),'InitialMagnification',100);
%     fig = figure; montage({outputImageLabel,outputImage}, 'Size', [1 2]);
%     set(fig,'Units','Inches');
%     pos = get(fig,'Position');
%     set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(fig,'D:\TCC\TCC\figs\im10_comparison.pdf','-dpdf','-r0');

    
%     fig = figure; montage({originalRgbImage,outputImageLabel,outputImage}, 'Size', [1 3]);
    
    figure; montage({rgbImage,outputImage});
end
classesTemp = classesTempJoin;
Kc = KcJoin;
centroids = centroidsJoin;
Uc = UcJoin;

% testSegment;

% end
% classes = classesTemp;
% colorSuperpixel;

% fig = figure; plot(2:30, jotas, 'linewidth', 2);
% xlabel('Nº of Clusters', 'fontsize', 15); xlim([2 30]);
% ylabel('Objective Function', 'fontsize', 15);
% title('Objective Function Values x Number of Clusters', 'fontsize', 20);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',15);

dbg = 1;
