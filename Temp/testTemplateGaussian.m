close all
allTemplateResults = [{}];
load R.mat;
R = R_copy;

for imgCnt = 1:10

%     if imgCnt == 3 || imgCnt == 10
%         dbg = 1;
%     else
%         continue;
%     end
    
load(['matriz_corre_img' num2str(imgCnt) '.mat']); 

%% resultado real
% load R.mat;
results = zeros(size(posMatrix));
% currResult = 2;
for k = 1:size(posMatrix, 1)
    for j = 1:size(posMatrix, 2)
        [lat_res, lon_res] = pix2latlon(R, posMatrix{k, j}(1), posMatrix{k, j}(2));
        results(k,j) = m_idist(lon_res, lat_res, lon(imgCnt), lat(imgCnt));
    end
end
[droneX, droneY] = find(results == min(min(results)), 1);
drone = posMatrix{droneX, droneY};
% resRes = 

% expandedResults
expandedPos = [linspace(cropSize(1), cropSize(3), 100)' linspace(cropSize(2), cropSize(4), 100)'];
expandedResults = zeros(100);
for k = 1:100
    for j = 1:100
        [lat_res, lon_res] = pix2latlon(R, expandedPos(k, 1), expandedPos(j, 2));
        expandedResults(k,j) = m_idist(lon_res, lat_res, lon(imgCnt), lat(imgCnt));
    end
end

% mapa completo
% plotMask(maskGeo, classesGeo, resultsMatrix{droneX, droneY}.rec, [], [], LmaskGeo, cropSize, drone);

% grid resultados
% figure; imagesc(results); colorbar;
% title('Real Density Matrix');
    
%% Resultado algoritmo

% close all;
% matriz de resultados
% load matriz_corre_img5.mat
nP = size(resultsMatrix, 1);

% 1 vetor só
allOut = [];
gridCnt = 1;
for k = 1:size(resultsMatrix, 1)
    for j = 1:size(resultsMatrix, 2)
        
        if isfield(resultsMatrix{k, j},'out')
%             if imgCnt == 1
                allOut = [allOut resultsMatrix{k, j}.out];
%             else
%                 allOut = [allOut resultsMatrix{k, j}.out(107)];
%             end
        else
            allOut = [allOut zeros(107, 1)];
        end
        
%         plotMask(resultsMatrix{k, j}.sub.mask, resultsMatrix{k, j}.sub.class,...
%         resultsMatrix{k, j}.sub.rec, [], [], resultsMatrix{k, j}.sub.Lmask);
%         title(['Figure X = ' num2str(j) ', Y = ' num2str(k) ', Cnt = ' num2str(gridCnt)]);
%         figure; montage({BW_Template, resultsMatrix{k, j}.sub.targ});
%         title(['Figure X = ' num2str(j) ', Y = ' num2str(k) ', Cnt = ' num2str(gridCnt)]);
%         gridCnt = gridCnt+1;
        
    end
end
% [~, idxMax] = max(allOut, [], 2);
% goodResultMax = sum(idxMax == [23 24 25 26]); % imagem 1, 24/25 é o melhor
% goodResultMax = sum(idxMax == [36 37 38 39 45 46]); % imagem 2, 39 é o melhor
% goodResultMax = sum(idxMax == [24 25 26 27 28 30 31 32 33 34 35]); % imagem 5, X é o melhor

% [~, idxMin] = min(allOut, [], 2);
% goodResultMin = sum(idxMin == [23 24 25 26]); % imagem 1, 25 é o melhor
% goodResultMin = sum(idxMin == [36 37 38 39 45 46]); % imagem 2, 39 é o melhor
% goodResultMin = sum(idxMin == [24 25 26 27 28 30 31 32 33 34 35]); % imagem 5, X é o melhor

bestVisPos = [55 66 64 45 64 35 47 56 57 61];
% 
% % jeito novo
idxCorreto = bestVisPos(imgCnt);
yc = ceil(idxCorreto/nP);
xc = mod(idxCorreto, nP);
% yc = 1; xc = 1;
visPos = posMatrix{yc, xc};

% matriz resultados
simMatrix = zeros(nP, nP);
metrica = 1;
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

[yMinRes, xMinRes] = find(expandedResults == min(expandedResults(:)), 1);
bestResult = [expandedPos(yMinRes, 1) expandedPos(xMinRes, 2)];
errY = abs(yCoords-bestResult(1));
yRes = find(errY == min(errY), 1);
errX = abs(xCoords-bestResult(2));
xRes = find(errX == min(errX), 1);

% plot distribution
z = z./max(z(:));
% figure;
% surf(x11,x22,z);
% figure;
% surf(yCoords,xCoords,z);
% title('UAV Position Distribution', 'fontsize', 20);
% xlabel('X Pixels', 'fontsize',15);
% ylabel('Y Pixels', 'fontsize', 15);
% zlabel('Probability', 'fontsize', 15);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',15);
% hold on; plot3(bestResult(1), bestResult(2), 1, 'r*', 'Linewidth', 2);


zCut = z;
zCut(zCut > 0.3) = 1;
zCut(zCut <= 0.3) = 0;
% figure;
% surf(x11,x22,zCut);
% figure;
% surf(yCoords,xCoords,zCut);
% imshow(flip(zCut));
% title('Potential UAV Area', 'fontsize', 20);
% xlabel('X Pixels', 'fontsize',15);
% ylabel('Y Pixels', 'fontsize', 15);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',15);
% hold on; plot3(bestResult(1), bestResult(2), 1, 'r*', 'Linewidth', 2);
% hold on; plot(yRes, xRes, 'r*', 'Linewidth', 10);

zRes = z(xRes, yRes);
zResCut = zCut(xRes, yRes);

cutAreaY = X1(logical(zCut));
cutAreaX = X2(logical(zCut));
resArea = [min(cutAreaY) min(cutAreaX) max(cutAreaY) max(cutAreaX)];

if imgCnt < 11
    errY = abs(yCoords-visPos(1));
    yRes = find(errY == min(errY), 1);
    errX = abs(xCoords-visPos(2));
    xRes = find(errX == min(errX), 1);

    zVis = z(xRes, yRes);
    zVisCut = zCut(xRes, yRes);
else
    zVis = 1e5;
    zVisCut = 0;
end

dbg = 1;

% grid resultados
% figure; imagesc(simMatrix); colorbar;
% title('Similarity Matrix');

resCentro = mu;
resProb = zRes;
resIn = zResCut;
resArea = sum(zCut(:)/numel(zCut));
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

%% results
% checkTemplateGridResults;

% Calculo da latitude e longitude com base nos pixeis da img.
[lat_res, lon_res] = pix2latlon(R, resCentro(1), resCentro(2));
resDist = m_idist(lon_res, lat_res, lon(imgCnt), lat(imgCnt));
if imgCnt < 11
    [lat_vis, lon_vis] = pix2latlon(R, visPos(1), visPos(2));
    visDist = m_idist(lon_res, lat_res, lon_vis, lat_vis);
else
    visDist = 1e5;
end

allTemplateResults = [allTemplateResults; {resProb, resIn, resDist, resArea, visProb, visIn, visDist}];


end

allTemplateResults
save templateResults_11.mat allTemplateResults
