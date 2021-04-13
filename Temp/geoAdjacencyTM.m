function [resCentro, resProb, resIn, resArea, resDist, visDist, visProb, visIn] = geoAdjacencyTM(geoAdjacencies, adjacencies, classes, Lmask, mask, maskGeo, cropSize, utilCropSize, parameters, filter, geo_img, util_mask, R, lon, lat, imgCnt, rotPlotImage, rotSegmentation)

%geoAdjacencies, ftGeoOwn, ftGeoAdj
run_template_matching = 1;
if run_template_matching
tic;
fprintf('\nRunning Correlation Matrix...\n');

% Extracting superpixels centroids and adjacency matrix
% tic;
sy = size(mask,1);
sx = size(mask,2);
removeSp = [];
sizesSp = zeros(1, length(adjacencies));
for labelVal = 1:length(adjacencies)
    [y, x] = find(mask == labelVal);
    if isempty(x)
        removeSp = [removeSp; labelVal];
        continue;
    end
    sizesSp(labelVal) = max(x)-min(x);
    adjacencies{labelVal}.centX = round(mean(x));
    adjacencies{labelVal}.centY = round(mean(y));
    adjacencies{labelVal}.cla = classes(labelVal);
    nX = min(x+1, sx);
    nY = min(y+1, sy);
    pX = max(x-1, 1);
    pY = max(y-1, 1);
    adjacents = zeros(length(x), 4);
    adjacents(:, 1) = mask(sub2ind(size(mask), y, nX));
    adjacents(:, 2) = mask(sub2ind(size(mask), nY, x));
    adjacents(:, 3) = mask(sub2ind(size(mask), y, pX));
    adjacents(:, 4) = mask(sub2ind(size(mask), pY, x));
    uniAdjs = unique(adjacents);
    uniAdjs = uniAdjs(uniAdjs > labelVal);
    adjacencies{labelVal}.adj = uniAdjs;
end

% removing any undesired superpixels
if ~isempty(removeSp)
    for labelVal = length(removeSp):-1:1
        adjacencies(removeSp(labelVal)) = [];
    end
end
% toc
dgb = 1;

% pegando retangulo do meio
x0 = max(utilCropSize(2), 1);
y0 = max(utilCropSize(1), 1);
x1 = min(utilCropSize(4), size(mask,2));
y1 = min(utilCropSize(3), size(mask,1));

% filtrando a parte útil
sy = size(Lmask,1);
sx = size(Lmask,2);
spCent = zeros(length(adjacencies), 2);
for spix = 1:length(adjacencies)
    py = adjacencies{spix}.centY;
    px = adjacencies{spix}.centX;
    if adjacencies{spix}.centY < utilCropSize(1) || adjacencies{spix}.centY > utilCropSize(3) ||...
       adjacencies{spix}.centX < utilCropSize(2) || adjacencies{spix}.centX > utilCropSize(4)
        px = Inf; py = Inf;
    end
    spCent(spix, :) = [py px];
end

rec = [y0 x0 y1 x1];
plotMask(mask, classes, rec);
figure; imagesc(util_mask);
pause(0.001);

BW_Template = edge(util_mask, 'canny', [], sqrt(2));

% finding class adjacency
% adjacencies,ftOwn,ftAdj
syG = size(maskGeo, 1);
sxG = size(maskGeo, 2);
tYG = cropSize(3)-cropSize(1);
tXG = cropSize(4)-cropSize(2);
nP = 10;
dy = round(tYG/(2*nP));
dx = round(tXG/(2*nP));
% nSP = length(adjacencies);
spCent = zeros(length(geoAdjacencies), 2);
corrMatrix = zeros(nP, nP);
corrMatrixEuc = zeros(nP, nP);
posMatrix = cell(nP);
resultsMatrix = cell(nP);
sizesSpGeo = zeros(1, length(geoAdjacencies));
for spix = 1:length(geoAdjacencies)
    py = geoAdjacencies{spix}.centY;
    px = geoAdjacencies{spix}.centX;
    sizesSpGeo(spix) = geoAdjacencies{spix}.sizeX;
    if geoAdjacencies{spix}.centY < cropSize(1) || geoAdjacencies{spix}.centY > cropSize(3) ||...
       geoAdjacencies{spix}.centX < cropSize(2) || geoAdjacencies{spix}.centX > cropSize(4)
        px = Inf; py = Inf;
    end
    spCent(spix, :) = [py px];
end

% Histogram equalization
labImage = rgb2lab(geo_img);
L = labImage(:,:,1)/100;
L = adapthisteq(L);
labImage(:,:,1) = L*100;
geo_img = lab2rgb(labImage);

pY = cropSize(1)+dy;
for ponY = 1:nP
    fprintf('\nLinha %d\n', ponY);
    pX = cropSize(2)+dx;
    for ponX = 1:nP
        fprintf('Coluna %d\n', ponX);
        
        posMatrix{ponY, ponX} = [pY pX];
        
        % pegando sub imagem da georreferenciada
        y0 = round(pY-sy/2); y1 = round(pY+sy/2); yc = 0;
        if y0 < 1; yc = -(1-y0); y0 = 1; y1 = sy; elseif y1 > syG; yc = y1-syG; y0 = syG-sy; y1 = syG; end
        x0 = round(pX-sx/2); x1 = round(pX+sx/2); xc = 0;
        if x0 < 1; xc = -(1-x0); x0 = 1; x1 = sx; elseif x1 > sxG; xc = x1-sxG; x0 = sxG-sx; x1 = sxG; end
        subImg = geo_img(y0:y1, x0:x1, :);
        [newGeoAdjacencies, ~, ~, maskGeoCrop, classesGeoCrop, idxGeoCrop] = getAdjacencies(subImg, parameters, filter);
        title(['Figure X = ' num2str(ponX) ', Y = ' num2str(ponY)]);
        
%         [outputImage1] = evalFunction(classes, length(unique(classes)), maskIdx, plotImg, length(classes));
%         [outputImage2] = evalFunction(classesGeoCrop, length(unique(classesGeoCrop)), idxGeoCrop, subImg, length(classesGeoCrop));
%         fig = figure; montage({subImg, outputImage2, rotPlotImage, rotSegmentation}, 'Size', [2 2]);
        
        spCent = zeros(length(newGeoAdjacencies), 2);
        sizesSpGeo = zeros(1, length(newGeoAdjacencies));
        oldGeoAdjacencies = newGeoAdjacencies;
        for adj = 1:length(newGeoAdjacencies)
            newGeoAdjacencies{adj}.centX = newGeoAdjacencies{adj}.centX + x0-1;
            newGeoAdjacencies{adj}.centY = newGeoAdjacencies{adj}.centY + y0-1;
            spCent(adj, :) = [newGeoAdjacencies{adj}.centY newGeoAdjacencies{adj}.centX];
            sizesSpGeo(adj) = newGeoAdjacencies{adj}.sizeX;
            newGeoAdjacencies{adj}.cla = classesGeoCrop(adj);
        end
        geoAdjacencies = newGeoAdjacencies;
        
        % colocando grid retangular na área útil
        tY = utilCropSize(3)-utilCropSize(1);
        tX = utilCropSize(4)-utilCropSize(2);
        currSpCent = spCent;
        for spix = 1:length(spCent)
            if currSpCent(spix, 1) ~= Inf
                if geoAdjacencies{spix}.centY < pY-(tY/2) || geoAdjacencies{spix}.centY > pY+(tY/2) ||...
                   geoAdjacencies{spix}.centX < pX-(tX/2) || geoAdjacencies{spix}.centX > pX+(tX/2)
                    currSpCent(spix, :) = [Inf Inf];
                end
            end
        end
        
        recGeo = [max(pY-(tY/2), 1) max(pX-(tX/2), 1) min(pY+(tY/2), syG) min(pX+(tX/2), sxG)];
                
        resultsMatrix{ponY, ponX}.rec = recGeo;
        
%         plotMask(maskGeo, classesGeo, recGeo, [], [], LmaskGeo, cropSize);
        LmaskGeoCrop = zeros(size(maskGeoCrop));
        for mx = 1:length(classesGeoCrop)
            LmaskGeoCrop(maskGeoCrop == mx) = classesGeoCrop(mx);
        end
        
        y0 = (sy/2)+yc-(tY/2); y1 = (sy/2)+yc+(tY/2);
        if y0 < 1; y0 = 1; y1 = tY+1; elseif y1 > size(subImg, 1); y0 = size(subImg, 1)-tY-1; y1 = size(subImg, 1); end
        x0 = (sx/2)+xc-(tX/2); x1 = (sx/2)+xc+(tX/2);
        if x0 < 1; x0 = 1; x1 = tX+1; elseif x1 > size(subImg, 2); x0 = size(subImg, 2)-tX-1; x1 = size(subImg, 2); end
%         recGeoCrop = [max((sy/2)+yc-(tY/2)-inc*sy, 1) max((sx/2)+xc-(tX/2)-inc*sx, 1) min((sy/2)+yc+(tY/2)+inc*sy, sy) min((sx/2)+xc+(tX/2)+inc*sx, sx)];
        
        % tendo certeza da dimensão igual
        if y1-y0+1 ~= size(util_mask, 1)
            shiftY = (y1-y0+1 - size(util_mask, 1));
            if y0-shiftY < 1; y1 = y1-shiftY; else; y0 = y0-shiftY; end
        end
        if x1-x0+1 ~= size(util_mask, 2)
            shiftX = (x1-x0+1 - size(util_mask, 2));
            if x0-shiftX < 1; x1 = x1-shiftX; else; x0 = x0-shiftX; end
        end
        
        recGeoCrop = [y0 x0 y1 x1];
%         plotMask(maskGeoCrop, classesGeoCrop, recGeoCrop, [], [], LmaskGeoCrop);    
%         title(['Figure X = ' num2str(ponX) ', Y = ' num2str(ponY)]);
        
        %subImage.img = subImg;
        subImage.rec = recGeoCrop;
        subImage.mask = maskGeoCrop;
        subImage.class = classesGeoCrop;
        subImage.Lmask = LmaskGeoCrop;
        
%         geo_mask = maskGeoCrop(recGeoCrop(1):recGeoCrop(3), recGeoCrop(2):recGeoCrop(4));
        util_geo_mask = LmaskGeoCrop(recGeoCrop(1):recGeoCrop(3), recGeoCrop(2):recGeoCrop(4));
        BW_Target = edge(util_geo_mask, 'canny', [], sqrt(2));
%         figure; montage({BW_Template, BW_Target});
%         title(['Figure X = ' num2str(ponX) ', Y = ' num2str(ponY)]);
        
        subImage.targ = BW_Target;
        resultsMatrix{ponY, ponX}.sub = subImage;
        
        BW_Template_mean=BW_Template-mean(mean(BW_Template));
        BW_Target_mean=BW_Target-mean(mean(BW_Target));

        try
            corrMat = xcorr2(BW_Target_mean,BW_Template_mean);
%             resultsMatrix{ponY, ponX}.out = [simbin(BW_Template,BW_Target); max(corrMat(:))];
            resultsMatrix{ponY, ponX}.out = max(corrMat(:));
        catch
            dbg = 1;
        end
        
%         figure; mesh(corrMat);
%         resultsMatrix{ponY, ponX}.corr = max(corrMat(:));
        
        save(['matriz_corre_img' num2str(imgCnt) '.mat'], 'posMatrix', 'resultsMatrix', 'pX', 'pY', 'ponY', 'ponX', 'BW_Template', 'cropSize');
        
        pX = pX + 2*dx;
    end
    pY = pY + 2*dy;
end
fprintf('\Tempo total = %f\n\n', toc);
save(['matriz_corre_img' num2str(imgCnt) '.mat'], 'posMatrix', 'resultsMatrix', 'pX', 'pY', 'ponY', 'ponX', 'BW_Template', 'cropSize');
dbg = 1;
else
   load(['matriz_corre_img' num2str(imgCnt) '.mat']); 
end

%% resultado real
% load R.mat;
results = zeros(size(posMatrix));
% currResult = 2;
for k = 1:size(posMatrix, 1)
    for j = 1:size(posMatrix, 2)
        [lat_res, lon_res] = pix2latlon(R, posMatrix{k, j}(1), posMatrix{k, j}(2));
        results(k,j) = m_idist(lon_res, lat_res, lon, lat);
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
        expandedResults(k,j) = m_idist(lon_res, lat_res, lon, lat);
    end
end

% mapa completo
% plotMask(maskGeo, classesGeo, resultsMatrix{droneX, droneY}.rec, [], [], LmaskGeo, cropSize, drone);

% grid resultados
% fig = figure; imagesc(results); colorbar;
% fig.CurrentAxes.XTick = 1:10;
% title('Real Distance Matrix', 'fontsize', 20);
    
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
            allOut = [allOut resultsMatrix{k, j}.out];
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

% plot distribution
z = z./max(z(:));
% figure;
% surf(yCoords,xCoords,z);
% title('UAV Position Distribution', 'fontsize', 20);
% xlabel('X Pixels', 'fontsize',15);
% ylabel('Y Pixels', 'fontsize', 15);
% zlabel('Probability', 'fontsize', 15);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',15);

zCut = z;
zCut(zCut > 0.3) = 1;
zCut(zCut <= 0.3) = 0;
% figure;
% surf(yCoords,xCoords,zCut);
% title('Potential UAV Area', 'fontsize', 20);
% xlabel('X Pixels', 'fontsize',15);
% ylabel('Y Pixels', 'fontsize', 15);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',15);

[yMinRes, xMinRes] = find(expandedResults == min(expandedResults(:)), 1);
bestResult = [expandedPos(yMinRes, 1) expandedPos(xMinRes, 2)];
errY = abs(yCoords-bestResult(1));
yRes = find(errY == min(errY), 1);
errX = abs(xCoords-bestResult(2));
xRes = find(errX == min(errX), 1);

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
% fig = figure; imagesc(simMatrix); colorbar;
% fig.CurrentAxes.XTick = 1:10;
% title('Correlation Matrix', 'fontsize', 20);
% 
% xlabel('Grid Point', 'fontsize',15);
% ylabel('Grid Point', 'fontsize', 15);
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',15);



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

%% results
% checkTemplateGridResults;

% Calculo da latitude e longitude com base nos pixeis da img.
[lat_res, lon_res] = pix2latlon(R, resCentro(1), resCentro(2));
resDist = m_idist(lon_res, lat_res, lon, lat);
if imgCnt < 11
    [lat_vis, lon_vis] = pix2latlon(R, visPos(1), visPos(2));
    visDist = m_idist(lon_res, lat_res, lon_vis, lat_vis);
else
    visDist = 1e5;
end

end