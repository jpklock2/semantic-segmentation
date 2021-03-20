function [corrMatrix] = geoAdjacency(ftAdj, ftOwn, geoAdjacencies, adjacencies, ftGeoOwn, ftGeoAdj, classes, classesGeo, Lmask, LmaskGeo, mask, maskGeo, cropSize, utilCropSize, rgbImage, parameters, filter, geo_img)

%geoAdjacencies, ftGeoOwn, ftGeoAdj
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

% filtrando a parte útil
sy = size(Lmask,1);
sx = size(Lmask,2);
inc = 0.1;
spCent = zeros(length(adjacencies), 2);
for spix = 1:length(adjacencies)
    py = adjacencies{spix}.centY;
    px = adjacencies{spix}.centX;
    if adjacencies{spix}.centY < utilCropSize(1)-inc*sy || adjacencies{spix}.centY > utilCropSize(3)+inc*sy ||...
       adjacencies{spix}.centX < utilCropSize(2)-inc*sx || adjacencies{spix}.centX > utilCropSize(4)+inc*sx
        px = Inf; py = Inf;
    end
    spCent(spix, :) = [py px];
end

% selecionando distancias no kernel quadrado
tY = utilCropSize(3)-utilCropSize(1);
tX = utilCropSize(4)-utilCropSize(2);
pY = round(tY/2);
pX = round(tX/2);
[selectedDist, selectedIdx] = sort(pdist2(spCent, [pY pX], 'euclidean'));
nSP = length(selectedDist(selectedDist < Inf));
selectedAdj = adjacencies(selectedIdx(1:nSP));
selectedCent = spCent(selectedIdx(1:nSP), :);
selectedIdx = selectedIdx(1:nSP);
selectedSizes = sizesSp(selectedIdx(1:nSP));

% colocando na ordem
selectedAdjSort = [{}];
selectedCentSort = [];
selectedIdxSort = [];
stopFlag = 0;
swipe = median(selectedSizes);
while stopFlag == 0
    selectedX = selectedCent(selectedCent(:, 2) < swipe, :);
    selectedAdjX = selectedAdj(selectedCent(:, 2) < swipe);
    selectedI = selectedIdx(selectedCent(:, 2) < swipe);
    [~, sortedY] = sort(selectedX(:, 1));
    selectedCentSort = [selectedCentSort; selectedX(sortedY, :)];
    selectedAdjSort = [selectedAdjSort selectedAdjX(sortedY)];
    selectedIdxSort = [selectedIdxSort; selectedI(sortedY)];
    selectedAdj(selectedCent(:, 2) < swipe) = [];
    selectedIdx(selectedCent(:, 2) < swipe) = [];
    selectedCent(selectedCent(:, 2) < swipe, :) = [];
    if swipe == cropSize(4); stopFlag = 1; end
    swipe = swipe + median(selectedSizes);
    if swipe > cropSize(4); swipe = cropSize(4); end
end

% corrigindo adjacentes
for sp = 1:length(selectedAdjSort)
    remove = [];
    for adj = 1:length(selectedAdjSort{sp}.adj)
        newAdj = find(selectedIdxSort == selectedAdjSort{sp}.adj(adj), 1);
        if isempty(newAdj)
            remove = [remove; adj];
        else
            selectedAdjSort{sp}.adj(adj) = newAdj;
        end
    end
    selectedAdjSort{sp}.adj(remove) = [];
end

rec = [max(utilCropSize(1)-inc*sy, 1) max(utilCropSize(2)-inc*sx, 1) min(utilCropSize(3)+inc*sy, size(mask,1)) min(utilCropSize(4)+inc*sx, size(mask,2))];
plotMask(mask, classes, rec, selectedAdjSort, selectedIdxSort);
pause(0.001);

% montando matriz
ftMatrix = cell(nSP);
relations = sort(ftAdj(:, 1:2), 2);
adjCnt = 1;
colorsFts = [];
texturesFts = [];
for sp = 1:length(selectedAdjSort)
    spix = selectedIdxSort(sp);
    spClass = classes(spix);
    ftMatrix{sp, sp} = ftOwn(spix, :);
    colorsFts = [colorsFts; ftOwn(spix, :)];
    for adj = 1:length(selectedAdjSort{sp}.adj)
        currAdj = selectedAdjSort{sp}.adj(adj);
        currIdx = selectedIdxSort(currAdj);
        currClass = classes(currIdx);
%         if spClass ~= currClass
            x = find(relations(:, 1) == min(spix, currIdx) & relations(:, 2) == max(spix, currIdx));
            fts = extractFeaturesAdj(rgbImage, selectedAdjSort{sp}, selectedAdjSort{currAdj});
            ftMatrix{sp, currAdj} = fts;
            ftMatrix{currAdj, sp} = fts;
            texturesFts = [texturesFts; fts];
%             ftMatrix{sp, currAdj} = ftAdj(x, 3:end);
%             ftMatrix{currAdj, sp} = ftAdj(x, 3:end);
%         end
        adjCnt = adjCnt + 1;
    end
end

% normalizando matriz
cMean = mean(colorsFts);
cStd = std(colorsFts);
fMean = mean(texturesFts);
fStd = std(texturesFts);
for rw = 1:size(ftMatrix, 1)
    for cl = 1:size(ftMatrix, 2)
        if rw == cl
            ftMatrix{rw, cl} = (ftMatrix{rw, cl}-cMean)./(2.*(3.*cStd+1));
        else
            if ~isempty(ftMatrix{rw, cl})
                ftMatrix{rw, cl} = (ftMatrix{rw, cl}-fMean)./(2.*(3.*fStd+1));
            end
        end
    end
end

% finding class adjacency
% adjacencies,ftOwn,ftAdj
% ftMatrix = cell(length(adjacencies));
% relations = sort(ftAdj(:, 1:2), 2);
% adjCnt = 1;
% for spix = 1:length(adjacencies)
%     spClass = classes(spix);
%     ftMatrix{spix, spix} = ftOwn(spix, :);
%     for adj = 1:length(adjacencies{spix}.adj)
%         currAdj = adjacencies{spix}.adj(adj);
%         currClass = classes(currAdj);
%         if spClass ~= currClass
%             x = find(relations(:, 1) == min(spix, currAdj) & relations(:, 2) == max(spix, currAdj));
%             ftMatrix{spix, currAdj} = ftAdj(x, 3:end);
%             ftMatrix{currAdj, spix} = ftAdj(x, 3:end);
%         end
%         adjCnt = adjCnt + 1;
%         dbg = 1;
%     end
% end
% dbg = 1;

% finding class adjacency
% adjacencies,ftOwn,ftAdj
syG = size(maskGeo, 1);
sxG = size(maskGeo, 2);
tYG = cropSize(3)-cropSize(1);
tXG = cropSize(4)-cropSize(2);
incG = 0.2;
cropSize(1) = max(1, cropSize(1)-tYG*incG);
cropSize(2) = max(1, cropSize(2)-tXG*incG);
cropSize(3) = min(syG, cropSize(3)+tYG*incG);
cropSize(4) = min(sxG, cropSize(4)+tXG*incG);
nP = 3;
dy = round((tYG+tYG*2*incG)/(2*nP));
dx = round((tXG+tXG*2*incG)/(2*nP));
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

incG = 0.5;
pY = cropSize(1)+dy;
for ponY = 1:nP
    fprintf('Linha %d\n', ponY);
    pX = cropSize(2)+dx;
    for ponX = 1:nP
        
        posMatrix{ponY, ponX} = [pY pX];
        
        % pegando sub imagem da georreferenciada
        y0 = round(pY-sy/2); y1 = round(pY+sy/2); yc = 0;
        if y0 < 1; yc = -(1-y0); y0 = 1; y1 = sy; elseif y1 > syG; yc = y1-syG; y0 = syG-sy; y1 = syG; end
        x0 = round(pX-sx/2); x1 = round(pX+sx/2); xc = 0;
        if x0 < 1; xc = -(1-x0); x0 = 1; x1 = sx; elseif x1 > sxG; xc = x1-sxG; x0 = sxG-sx; y1 = sxG; end
        subImg = geo_img(y0:y1, x0:x1, :);
        [newGeoAdjacencies, newFtGeoOwn, newFtGeoAdj, maskGeoCrop, classesGeoCrop] = getAdjacencies(subImg, parameters, filter);
        spCent = zeros(length(newGeoAdjacencies), 2);
        geoAdjacenciesCrop = newGeoAdjacencies;
        sizesSpGeo = zeros(1, length(newGeoAdjacencies));
        for adj = 1:length(newGeoAdjacencies)
            newGeoAdjacencies{adj}.centX = newGeoAdjacencies{adj}.centX + x0;
            newGeoAdjacencies{adj}.centY = newGeoAdjacencies{adj}.centY + y0;
            spCent(adj, :) = [newGeoAdjacencies{adj}.centY newGeoAdjacencies{adj}.centX];
            sizesSpGeo(adj) = newGeoAdjacencies{adj}.sizeX;
        end
        geoAdjacencies = newGeoAdjacencies;
        ftGeoOwn = newFtGeoOwn;
        ftGeoAdj = newFtGeoAdj;
        
        % colocando grid retangular na área útil
        tY = utilCropSize(3)-utilCropSize(1);
        tX = utilCropSize(4)-utilCropSize(2);
        currSpCent = spCent;
        for spix = 1:length(spCent)
            if currSpCent(spix, 1) ~= Inf
                if geoAdjacencies{spix}.centY < pY-(tY/2)-inc*sy || geoAdjacencies{spix}.centY > pY+(tY/2)+inc*sy ||...
                   geoAdjacencies{spix}.centX < pX-(tX/2)-inc*sx || geoAdjacencies{spix}.centX > pX+(tX/2)+inc*sx
                    currSpCent(spix, :) = [Inf Inf];
                end
            end
        end
        
        % colocando grid retangular
%         currSpCent = spCent;
%         for spix = 1:length(spCent)
%             if currSpCent(spix, 1) ~= Inf
%                 if geoAdjacencies{spix}.centY < pY-(tY/2)-incG*sy || geoAdjacencies{spix}.centY > pY+(tY/2)+incG*sy ||...
%                    geoAdjacencies{spix}.centX < pX-(tX/2)-incG*sx || geoAdjacencies{spix}.centX > pX+(tX/2)+incG*sx
%                     currSpCent(spix, :) = [Inf Inf];
%                 end
%             end
%         end
        
        [~, selectedIdx] = sort(pdist2(currSpCent, [pY pX], 'euclidean'));
        selectedAdj = geoAdjacencies(selectedIdx(1:nSP));
        selectedAdjCrop = geoAdjacenciesCrop(selectedIdx(1:nSP));
        selectedCent = currSpCent(selectedIdx(1:nSP), :);
        selectedIdx = selectedIdx(1:nSP);
        selectedSizes = sizesSpGeo(selectedIdx(1:nSP));
        
        % colocando na ordem
        selectedAdjSort = [{}];
        selectedAdjSortCrop = [{}];
        selectedCentSort = [];
        selectedIdxSort = [];
        stopFlag = 0;
%         swipe = cropSize(2)+floor((tX+tX*2*incG)/25);
        swipe = median(selectedSizes);
        while stopFlag == 0
            selectedX = selectedCent(selectedCent(:, 2) < swipe, :);
            selectedAdjX = selectedAdj(selectedCent(:, 2) < swipe);
            selectedAdjXCrop = selectedAdjCrop(selectedCent(:, 2) < swipe);
            selectedI = selectedIdx(selectedCent(:, 2) < swipe);
            [~, sortedY] = sort(selectedX(:, 1));
            selectedCentSort = [selectedCentSort; selectedX(sortedY, :)];
            selectedAdjSort = [selectedAdjSort selectedAdjX(sortedY)];
            selectedAdjSortCrop = [selectedAdjSortCrop selectedAdjXCrop(sortedY)];
            selectedIdxSort = [selectedIdxSort; selectedI(sortedY)];
            selectedAdj(selectedCent(:, 2) < swipe) = [];
            selectedAdjCrop(selectedCent(:, 2) < swipe) = [];
            selectedIdx(selectedCent(:, 2) < swipe) = [];
            selectedCent(selectedCent(:, 2) < swipe, :) = [];
            if swipe == cropSize(4); stopFlag = 1; end
%             swipe = swipe + floor((tX+tX*2*incG)/25);
            swipe = swipe + median(selectedSizes);
            if swipe > cropSize(4); swipe = cropSize(4); end
        end
        
        % corrigindo adjacentes
        for sp = 1:length(selectedAdjSort)
            remove = [];
            for adj = 1:length(selectedAdjSort{sp}.adj)
                newAdj = find(selectedIdxSort == selectedAdjSort{sp}.adj(adj), 1);
                if isempty(newAdj)
                    remove = [remove; adj];
                else
                    selectedAdjSort{sp}.adj(adj) = newAdj;
                    selectedAdjSortCrop{sp}.adj(adj) = newAdj;
                end
            end
            selectedAdjSort{sp}.adj(remove) = [];
            selectedAdjSortCrop{sp}.adj(remove) = [];
            if sp == 30
                dbg = 1;
            end
        end
        
%         recGeo = [max(pY-(tY/2)-incG*sy, cropSize(1)) max(pX-(tX/2)-incG*sx, cropSize(2)) min(pY+(tY/2)+incG*sy, cropSize(3)) min(pX+(tX/2)+incG*sx, cropSize(4))];
        recGeo = [max(pY-(tY/2)-inc*sy, 1) max(pX-(tX/2)-inc*sx, 1) min(pY+(tY/2)+inc*sy, syG) min(pX+(tX/2)+inc*sx, sxG)];
                
        resultsMatrix{ponY, ponX}.rec = recGeo;
        resultsMatrix{ponY, ponX}.adj = selectedAdjSort;
        resultsMatrix{ponY, ponX}.idx = selectedIdxSort;
        
        plotMask(maskGeo, classesGeo, recGeo, selectedAdjSort, selectedIdxSort, LmaskGeo, cropSize);
        LmaskGeoCrop = zeros(size(maskGeoCrop));
        for mx = 1:length(classesGeoCrop)
            LmaskGeoCrop(maskGeoCrop == mx) = classesGeoCrop(mx);
        end
        recGeoCrop = [max((sy/2)+yc-(tY/2)-inc*sy, 1) max((sx/2)+xc-(tX/2)-inc*sx, 1) min((sy/2)+yc+(tY/2)+inc*sy, sy) min((sx/2)+xc+(tX/2)+inc*sx, sx)];
        
        subImage.img = subImg;
        subImage.rec = recGeoCrop;
        subImage.mask = maskGeoCrop;
        subImage.class = classesGeoCrop;
        subImage.adj = selectedAdjSortCrop;
        subImage.idx = selectedIdxSort;
        subImage.Lmask = LmaskGeoCrop;
        
        resultsMatrix{ponY, ponX}.sub = subImage;
        
        plotMask(maskGeoCrop, classesGeoCrop, recGeoCrop, selectedAdjSortCrop, selectedIdxSort, LmaskGeoCrop);    

        % montando matriz
        ftGeoMatrix = cell(nSP);
        relations = sort(ftGeoAdj(:, 1:2), 2);
        adjCnt = 1;
        for sp = 1:length(selectedAdjSort)
            spix = selectedIdxSort(sp);
            spClass = classesGeo(spix);
            ftGeoMatrix{sp, sp} = ftGeoOwn(spix, :);
            ftGeoMatrix{sp, sp} = (ftGeoMatrix{sp, sp}-cMean)./(2.*(3.*cStd+1));
            for adj = 1:length(selectedAdjSort{sp}.adj)
                currAdj = selectedAdjSort{sp}.adj(adj);
                currIdx = selectedIdxSort(currAdj);
                currClass = classesGeo(currIdx);
%                 if spClass ~= currClass
                    x = find(relations(:, 1) == min(spix, currIdx) & relations(:, 2) == max(spix, currIdx));
                    ftGeoMatrix{sp, currAdj} = ftGeoAdj(x, 3:end);
                    ftGeoMatrix{currAdj, sp} = ftGeoAdj(x, 3:end);
                    ftGeoMatrix{sp, currAdj} = (ftGeoMatrix{sp, currAdj}-fMean)./(2.*(3.*fStd+1));
                    ftGeoMatrix{currAdj, sp} = (ftGeoMatrix{currAdj, sp}-fMean)./(2.*(3.*fStd+1));
%                 end
                adjCnt = adjCnt + 1;
            end
        end
        
        % calculando distancia
        for p1 = 1:size(ftMatrix, 1)
            for p2 = p1:size(ftMatrix, 2)
%                 disp([ponY ponX p1 p2]);
                z = 0;
                zg = 0;
                if p1 ~= p2; currFt = zeros(1, 270); currGeoFt = zeros(1, 270); end
                if ~isempty(ftMatrix{p1,p2}); currFt = ftMatrix{p1,p2}; z = 1; end
                if ~isempty(ftGeoMatrix{p1,p2}); currGeoFt = ftGeoMatrix{p1,p2}; zg = 1; end
                if (z == 0 &&  zg == 1) || (z == 1 &&  zg == 0)
                    corrMatrix(ponY, ponX) = corrMatrix(ponY, ponX) + 1;
                elseif (z == 1 && zg == 1)
                    corrMatrix(ponY, ponX) = corrMatrix(ponY, ponX) + pdist2(currFt, currGeoFt, 'cosine');
                end
                try
                    corrMatrixEuc(ponY, ponX) = corrMatrixEuc(ponY, ponX) + pdist2(currFt, currGeoFt, 'euclidean');
                catch
                    dbg = 1;
                end
            end
        end
        
        pX = pX + 2*dx;
    end
    pY = pY + 2*dy;
end
fprintf('\Tempo total = %f\n\n', toc);
save matriz_corre_22_01.mat corrMatrix corrMatrixEuc posMatrix resultsMatrix
dbg = 1;

% matriz de resultados
load matriz_corre_22_01.mat
load R.mat;
results = zeros(size(posMatrix));
for k = 1:size(posMatrix, 1)
    for j = 1:size(posMatrix, 2)
        [lat_srp, lon_srp] = pix2latlon(R_copy, posMatrix{k, j}(1), posMatrix{k, j}(2));
        results(k,j) = m_idist(lon_srp, lat_srp, lon(1,1), lat(1,1));
    end
end

[droneX, droneY] = find(results == min(min(results)), 1);
drone = posMatrix{droneX, droneY};
% ideal result
plotMask(maskGeo, classesGeo, resultsMatrix{droneX, droneY}.rec, resultsMatrix{droneX, droneY}.adj, resultsMatrix{droneX, droneY}.idx, LmaskGeo, cropSize, drone);

% jeito novo
[yc, xc] = find(corrMatrix == min(min(corrMatrix)));
fprintf('\nResult for Cosine Distance %.2f\n', results(yc(1), xc(1)));
[ye, xe] = find(corrMatrixEuc == min(min(corrMatrixEuc)));
fprintf('Result for Euclidean Distance %.2f\n', results(ye(1), xe(1)));

droneX = yc;
droneY = xc;
plotMask(maskGeo, classesGeo, resultsMatrix{droneX, droneY}.rec, resultsMatrix{droneX, droneY}.adj, resultsMatrix{droneX, droneY}.idx, LmaskGeo, cropSize, drone);

% resultsCorr = zeros(size(corrMatrix));
% resultsCorrEuc = zeros(size(corrMatrix));
% for k = 1:size(corrMatrix, 1)
%     for j = 1:size(corrMatrix, 2)
%         [lat_srp, lon_srp] = pix2latlon(R_copy, corrMatrix(k, j)(1), corrMatrix{k, j}(2));
%         resultsCorr(k,j) = m_idist(lon_srp, lat_srp, lon(1,1), lat(1,1));
%         [lat_srp, lon_srp] = pix2latlon(R_copy, resultsCorrEuc{k, j}(1), resultsCorrEuc{k, j}(2));
%         resultsCorrEuc(k,j) = m_idist(lon_srp, lat_srp, lon(1,1), lat(1,1));
%     end
% end

% distM = pdist2(textureCentroids, textureCentroids, 'correlation');
figure; imagesc(results); colorbar;
title('Real Density Matrix');

figure; imagesc(corrMatrix);
title('Correlation Cosine Density Matrix');

figure; imagesc(corrMatrixEuc);
title('Correlation Euclidean Density Matrix');

plotMask(resultsMatrix{5, 5}.sub.mask, resultsMatrix{5, 5}.sub.class,...
        resultsMatrix{5, 5}.sub.rec, resultsMatrix{5, 5}.sub.adj,...
        resultsMatrix{5, 5}.sub.idx, resultsMatrix{5, 5}.sub.Lmask);

        subImage.img = subImg;
        subImage.rec = recGeoCrop;
        subImage.mask = maskGeoCrop;
        subImage.class = classesGeoCrop;
        subImage.adj = selectedAdjSortCrop;
        subImage.idx = selectedIdxSort;
        subImage.Lmask = LmaskGeoCrop;
        
        resultsMatrix{ponY, ponX}.sub = subImage;

% M=size(corrMatrix);
% M2=resultsCorr;
% mat_string=reshape(arrayfun(@(x) strtrim(cellstr(sprintf('%2.2f',x))),M2),M);
% X = 1:size(corrMatrix, 1);
% Y = 1:size(corrMatrix, 2);
% for i = 1:M(1)
%     for j = 1:M(2)
%         if (corrMatrix(i,j) < 1)
%             text(X(j),Y(i),mat_string{i,j},'horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','white');
%         else
%             text(X(j),Y(i),mat_string{i,j},'horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','black');
%         end
%     end
% end


% jeito antigo
% [lat_srp, lon_srp] = pix2latlon(R_copy, posMatrix{12, 40}(1), posMatrix{12, 40}(2));
% m_idist(lon_srp, lat_srp, lon(i,1), lat(i,1))
