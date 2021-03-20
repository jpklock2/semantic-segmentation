close all;

Nexp = 100;
allUc = cell(Nexp, 1);
allJc = zeros(Nexp, 1);
allCentroids = cell(Nexp, 1);
for k = 1:Nexp
    [~, allUc{k}, JcTemp, allCentroids{k}] = MyFuzzyMeans_opt(pixelsOri, KF);
    allJc(k) = JcTemp(end);
end
medianJc = abs(allJc - median(allJc));
meanJc = abs(allJc - mean(allJc));
idxJcMedian = find(medianJc == min(medianJc), 1);
idxJcMean = find(meanJc == min(meanJc), 1);
idxJc = idxJcMedian;
idxJc = idxJcMean;
idxJc = find(allJc == min(allJc), 1);
Kc = KF;

figure;
images = [{}];
images = [images; rgbImage];
cophs = zeros(1, 3);
indName = [{}];

for indi = 1:3
    for met = 1:1

if indi == 1
    idxJc = idxJcMedian;
    indName = [indName; 'Median'];
elseif indi == 2
    idxJc = idxJcMean;
    indName = [indName; 'Mean'];
elseif indi == 3
    idxJc = find(allJc == min(allJc), 1);
    indName = [indName; 'Min'];
end

methods = {'centroid', 'median', 'ward', 'average', 'complete', 'single', 'weighted'};
metrics = {'euclidean', 'cityblock', 'minkowski', 'chebychev', 'cosine', 'correlation', 'spearman'};
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

if met == 1
    [maxMeth, maxMetr] = find(cDist == max(max(cDist)), 1);
% elseif met == 2
%     [maxMeth, maxMetr] = find(cClust == max(max(cClust)), 1);
% elseif met == 3
%     maxMeth = 4; maxMetr = 6;
end

cophs(indi) = cDist(maxMeth, maxMetr);
D = pdist(allCentroids{idxJc}, metrics{maxMetr});
tree = linkage(allCentroids{idxJc}, methods{maxMeth}, metrics{maxMetr});
leafOrder = optimalleaforder(tree, D);
% figure();
subplot(1,3,indi);
[H, T, OUTPERM] = dendrogram(tree, 'Reorder', leafOrder);
I = inconsistent(tree);
disp([indi met]);
disp(cDist(maxMeth, maxMetr));
disp(I(:, 4));

Uc = allUc{idxJc};
centroids = allCentroids{idxJc};
[~, idxUc] = sort(Uc, 2);
classesTemp = idxUc(:, end);
if plotsCompare
    thresh = 0.15;
    [classesTempJoin, KcJoin] = joinClasses(tree, Kc, classesTemp, thresh);
    [outputImage] = evalFunction(classesTempJoin, KcJoin, idx, rgbImage, N);
    images = [images; outputImage];
%     subplot(2,4,indi+5); imshow(outputImage);
%     figure; montage({rgbImage,outputImage});
%     title(['Index: ' indName ' - Cophenet: ' num2str(cDist(maxMeth, maxMetr))], 'FontSize', 16);
end

dbg = 1;
    end
end
suptitle('Dendrograms');
figure; montage(images);
set(gcf, 'units','normalized');
text(20, -20, 'Original Image', 'FontSize', 15);
text(1030, -20, ['Index: ' indName{1} ' - Cophenet: ' num2str(cophs(1))], 'FontSize', 15);
text(20, 1100, ['Index: ' indName{2} ' - Cophenet: ' num2str(cophs(2))], 'FontSize', 15);
text(1080, 1100, ['Index: ' indName{3} ' - Cophenet: ' num2str(cophs(3))], 'FontSize', 15);

% title(['Index: ' indName ' - Cophenet: ' num2str(cDist(maxMeth, maxMetr))], 'FontSize', 16);