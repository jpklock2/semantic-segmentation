load dendro.mat

medianJc = abs(allJc - median(allJc));
meanJc = abs(allJc - mean(allJc));
idxJcMedian = find(medianJc == min(medianJc), 1);
idxJcMean = find(meanJc == min(meanJc), 1);
idxJc = find(allJc == min(allJc), 1);

myMeth = 'weighted';
myMetric = 'jaccard';

tree = linkage(allCentroids{idxJc}, myMeth, myMetric);
D = pdist(allCentroids{idxJc}, myMetric);
leafOrder = optimalleaforder(tree, D);
figure();
[H, T, OUTPERM] = dendrogram(tree, 'Reorder', leafOrder);

title(['Method: ' myMeth '- Metric: ' myMetric]);

Uc = allUc{idxJc};
centroids = allCentroids{idxJc};
[~, idxUc] = sort(Uc, 2);
classesTemp = idxUc(:, end);
Kc = length(unique(classesTemp));
thresh = 0.2;
[classesTempJoin, KcJoin] = joinClasses(tree, Kc, classesTemp, thresh);
[outputImage] = evalFunction(classesTempJoin, KcJoin, idx, rgbImage, N);
figure; montage({rgbImage, outputImage});
