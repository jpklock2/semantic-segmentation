% close all;
% matriz de resultados
% load matriz_corre_img5.mat
nP = size(resultsMatrix, 1);

% 1 vetor só
allOut = [];
for k = 1:size(resultsMatrix, 1)
    for j = 1:size(resultsMatrix, 2)
        allOut = [allOut resultsMatrix{k, j}.out];
        
%         plotMask(resultsMatrix{k, j}.sub.mask, resultsMatrix{k, j}.sub.class,...
%         resultsMatrix{k, j}.sub.rec, [], [], resultsMatrix{k, j}.sub.Lmask);
%         title(['Figure X = ' num2str(j) ', Y = ' num2str(k)]);
%         figure; montage({BW_Template, resultsMatrix{k, j}.sub.targ});
%         title(['Figure X = ' num2str(j) ', Y = ' num2str(k)]);
        
    end
end
[~, idxMax] = max(allOut, [], 2);
% goodResultMax = sum(idxMax == [23 24 25 26]); % imagem 1, 25 é o melhor
% goodResultMax = sum(idxMax == [36 37 38 39 45 46]); % imagem 2, 39 é o melhor
goodResultMax = sum(idxMax == [24 25 26 27 28 30 31 32 33 34 35]); % imagem 5, X é o melhor

[~, idxMin] = min(allOut, [], 2);
% goodResultMin = sum(idxMin == [23 24 25 26]); % imagem 1, 25 é o melhor
% goodResultMin = sum(idxMin == [36 37 38 39 45 46]); % imagem 2, 39 é o melhor
goodResultMin = sum(idxMin == [24 25 26 27 28 30 31 32 33 34 35]); % imagem 5, X é o melhor

bestVisPos = [25 39];

% jeito novo
idxCorreto = 39;
yc = ceil(idxCorreto/nP);
xc = mod(idxCorreto, nP);
visPos = posMatrix{yc, xc};

% matriz resultados
simMatrix = zeros(nP, nP);
metrica = 107;
currResult = allOut(metrica, :);
cnt = 1;
for i = 1:nP
    for j = 1:nP
        simMatrix(i, j) = currResult(cnt);
        cnt = cnt+1;
    end
end

maxCorr = max(simMatrix(:));
threshCorr = 0.7*maxCorr;
coords = cell2mat(posMatrix(simMatrix >= threshCorr));

yCoords = linspace(cropSize(1), cropSize(3), 10000);
xCoords = linspace(cropSize(2), cropSize(4), 10000);
[X1,X2] = meshgrid(yCoords, xCoords);
X = [X1(:) X2(:)];

mu = mean(coords);
Sigma = cov(coords);
z = mvnpdf(X, mu, Sigma);
z = reshape(z,length(xCoords),length(yCoords));

% plot distribution
z = z./max(z(:));
% figure;
% surf(x11,x22,z);

zCut = z;
zCut(zCut > 0.3) = 1;
zCut(zCut <= 0.3) = 0;
% figure;
% surf(x11,x22,zCut);

[yMinRes, xMinRes] = find(expandedResults == min(expandedResults(:)), 1);
bestResult = [expandedPos(yMinRes, 1) expandedPos(xMinRes, 2)];
errY = abs(yCoords-bestResult(1));
yRes = find(errY == min(errY), 1);
errX = abs(xCoords-bestResult(2));
xRes = find(errX == min(errX), 1);

zRes = z(yRes, xRes);
zResCut = zCut(yRes, xRes);

cutAreaY = X1(logical(zCut));
cutAreaX = X2(logical(zCut));
resArea = [min(cutAreaY) min(cutAreaX) max(cutAreaY) max(cutAreaX)];

if imgCnt < 11
    errY = abs(yCoords-visPos(1));
    yRes = find(errY == min(errY), 1);
    errX = abs(xCoords-visPos(2));
    xRes = find(errX == min(errX), 1);

    zVis = z(yRes, xRes);
    zVisCut = zCut(yRes, xRes);
else
    zVis = 1e5;
    zVisCut = 0;
end

dbg = 1;

% grid resultados
figure; imagesc(simMatrix); colorbar;
title('Similarity Matrix');

resCentro = mu;
resProb = zRes;
resIn = zResCut;
visProb = zVis;
visIn = zVisCut;
% resArea = 
% idealArea = 1;

% mapa completo
% plotMask(resultsMatrix{yc, xc}.sub.mask, resultsMatrix{yc, xc}.sub.class,...
%         resultsMatrix{yc, xc}.sub.rec, [], [], resultsMatrix{yc, xc}.sub.Lmask);

% mascaras
% BW_Target = resultsMatrix{yc, xc}.sub.targ;
% figure; montage({BW_Template, BW_Target});
