function [corrMatrix] = geoAdjacencyGD(ftAdj, ftOwn, geoAdjacencies, adjacencies, ftGeoOwn, ftGeoAdj, classes, classesGeo, Lmask, LmaskGeo, mask, maskGeo, cropSize, utilCropSize, rgbImage, parameters, filter, geo_img, util_mask, util_image, centroidsGeo)

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
distribution = [];
for spix = 1:length(adjacencies)
    py = adjacencies{spix}.centY;
    px = adjacencies{spix}.centX;
    if adjacencies{spix}.centY < utilCropSize(1) || adjacencies{spix}.centY > utilCropSize(3) ||...
       adjacencies{spix}.centX < utilCropSize(2) || adjacencies{spix}.centX > utilCropSize(4)
        px = Inf; py = Inf;
    else
        distribution = [distribution; py-y0+1 px-x0+1 adjacencies{spix}.cla];
    end
    spCent(spix, :) = [py px];
end

rec = [y0 x0 y1 x1];

x11 = linspace(1, size(util_mask, 2), 100);
x22 = linspace(1, size(util_mask, 1), 100);
% [X1,X2] = meshgrid(x11, x22);
% X = [X1(:) X2(:)];
droneClasses = unique(distribution(:, 3));
% pdfs = [{}];
weigths = zeros(length(droneClasses), 1);
finalPdf = zeros(100, 100);
finalPdfNorm = zeros(100, 100);
pdfs = [{}];
for cla = 1:length(droneClasses)
    y = zeros(100, 100);
    class = droneClasses(cla);
    dist = distribution(distribution(:, 3) == class, 1:2);
    gm = getBestFit(dist, 0);
%     mu = mean(dist);
%     Sigma = cov(dist);
%     y = mvnpdf(X,mu,Sigma);
%     pdfs = [pdfs; y];
%     gm = gmdistribution(mu,Sigma);
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
    for j = 1:100
        y(:, j) = gmPDF(x22, x11(j)*ones(1, 100));
    end
    y = y';
%     normY = (y - min(y))./(max(y)-min(y));
    normY = (1/max(y(:)))*y;
    pdfs = [pdfs; normY];
    weigths(cla) = sum(sum(util_mask == class))/(size(util_mask, 1)*size(util_mask, 2));
    finalPdf = finalPdf + y*weigths(cla);
    finalPdfNorm = finalPdfNorm + normY*weigths(cla);
%     y = reshape(y,length(x22),length(x11));
%     figure;
%     surf(x22,x11,normY);
%     util_mask_temp = util_mask;
%     util_mask_temp(util_mask_temp ~= class) = 0;
%     figure; imagesc(util_mask_temp);
%     fsurf(gmPDF,[1 size(util_mask, 1) 1 size(util_mask, 2)]);
    dbg = 1;
end

% figure;
% surf(x22,x11,finalPdf);
% title('Final PDF');
% 
% figure;
% surf(x22,x11,finalPdfNorm);
% title('Normalized Final PDF');

plotMask(mask, classes, rec);
figure; imagesc(util_mask);
pause(0.001);

% BW_Template = edge(util_mask, 'canny', [], sqrt(2));

% finding class adjacency
% adjacencies,ftOwn,ftAdj
syG = size(maskGeo, 1);
sxG = size(maskGeo, 2);
tYG = cropSize(3)-cropSize(1);
tXG = cropSize(4)-cropSize(2);
incG = 0;
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

% distanceMatrix
% centDist = pdist2(centroidsGeo, centroidsGeo, 'euclidean');
centDist = pdist2(centroidsGeo, centroidsGeo, 'cosine');
centDist = 1-centDist./(max(centDist));

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
        [newGeoAdjacencies, ~, ~, maskGeoCrop, classesGeoCrop] = getAdjacencies(subImg, parameters, filter);
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
        geoDistribution = [];
        for spix = 1:length(spCent)
            if currSpCent(spix, 1) ~= Inf
                if geoAdjacencies{spix}.centY < pY-(tY/2) || geoAdjacencies{spix}.centY > pY+(tY/2) ||...
                   geoAdjacencies{spix}.centX < pX-(tX/2) || geoAdjacencies{spix}.centX > pX+(tX/2)
                    currSpCent(spix, :) = [Inf Inf];
                else
                    geoDistribution = [geoDistribution; oldGeoAdjacencies{spix}.centY oldGeoAdjacencies{spix}.centX geoAdjacencies{spix}.cla];
                end
            end
        end
        
        recGeo = [max(pY-(tY/2), 1) max(pX-(tX/2), 1) min(pY+(tY/2), syG) min(pX+(tX/2), sxG)];
                
        resultsMatrix{ponY, ponX}.rec = recGeo;
        
        plotMask(maskGeo, classesGeo, recGeo, [], [], LmaskGeo, cropSize);
        LmaskGeoCrop = zeros(size(maskGeoCrop));
        for mx = 1:length(classesGeoCrop)
            LmaskGeoCrop(maskGeoCrop == mx) = classesGeoCrop(mx);
        end
        
        y0 = (sy/2)+yc-(tY/2); y1 = (sy/2)+yc+(tY/2);
        if y0 < 1; y0 = 1; y1 = tY+1; elseif y1 > size(subImg, 1); y0 = size(subImg, 1)-tY-1; y1 = size(subImg, 1); end
        x0 = (sx/2)+xc-(tX/2); x1 = (sx/2)+xc+(tX/2);
        if x0 < 1; x0 = 1; x1 = tX+1; elseif x1 > size(subImg, 2); x0 = size(subImg, 2)-tX-1; x1 = size(subImg, 2); end
%         recGeoCrop = [max((sy/2)+yc-(tY/2)-inc*sy, 1) max((sx/2)+xc-(tX/2)-inc*sx, 1) min((sy/2)+yc+(tY/2)+inc*sy, sy) min((sx/2)+xc+(tX/2)+inc*sx, sx)];
        recGeoCrop = [y0 x0 y1 x1];
        
        for gdis = 1:size(geoDistribution, 1)
            geoDistribution(gdis, 1) = geoDistribution(gdis, 1) - recGeoCrop(1) + 1;
            geoDistribution(gdis, 2) = geoDistribution(gdis, 2) - recGeoCrop(2) + 1;
        end
        
        plotMask(maskGeoCrop, classesGeoCrop, recGeoCrop, [], [], LmaskGeoCrop);    
        
        %subImage.img = subImg;
        subImage.rec = recGeoCrop;
        subImage.mask = maskGeoCrop;
        subImage.class = classesGeoCrop;
        subImage.Lmask = LmaskGeoCrop;
        
        resultsMatrix{ponY, ponX}.sub = subImage;
        
        geo_mask = maskGeoCrop(recGeoCrop(1):recGeoCrop(3), recGeoCrop(2):recGeoCrop(4));
        util_geo_mask = LmaskGeoCrop(recGeoCrop(1):recGeoCrop(3), recGeoCrop(2):recGeoCrop(4));
        
        x11 = linspace(1, size(util_geo_mask, 2), 100);
        x22 = linspace(1, size(util_geo_mask, 1), 100);

        geoClasses = unique(geoDistribution(:, 3));

        geoPdfs = [{}];
        weigths = zeros(length(geoClasses), 1);
        geoFinalPdf = zeros(100, 100);
        geoFinalPdfNorm = zeros(100, 100);
        dtwResult = 0;
        dtwResult2 = 0;
        droneClassesTemp = droneClasses;
        
        for cla = 1:length(geoClasses)
            class = geoClasses(cla);
            dist = geoDistribution(geoDistribution(:, 3) == class, 1:2);
            if ~isempty(dist) && size(dist, 1) > 2
                y = zeros(100, 100);
                try
                    gm = getBestFit(dist, 0);
                    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(gm,[x0 y0]),x,y);
                catch
                    dbg2 = 1;
                end
                for j = 1:100
                    y(:, j) = gmPDF(x22, x11(j)*ones(1, 100));
                end
                y = y';

                normY = (1/max(y(:)))*y;
                geoPdfs = [geoPdfs; normY];
                weigths(cla) = sum(sum(util_geo_mask == class))/(size(util_geo_mask, 1)*size(util_geo_mask, 2));
                geoFinalPdf = geoFinalPdf + y*weigths(cla);
                geoFinalPdfNorm = geoFinalPdfNorm + normY*weigths(cla);

    %             figure;
    %             surf(x22,x11,normY);
    %             title(['Geo Image - Density for class ' num2str(class)]);
    %             util_geo_mask_temp = util_geo_mask;
    %             util_geo_mask_temp(util_geo_mask_temp ~= class) = 0;
    %             figure; imagesc(util_geo_mask_temp);
    %             title(['Geo Image - Mask for class ' num2str(class)]);

                for cla2 = 1:length(droneClasses)
                    dtwResult = dtwResult + dtw(pdfs{cla2}, normY, 100)*centDist(droneClasses(cla2), class);
                end

                if any(droneClasses == class)

    %                 figure;
    %                 surf(x22,x11,pdfs{find(droneClasses == class, 1)});
    %                 title(['Drone Image - Density for class ' num2str(class)]);
    %                 util_mask_temp = util_mask;
    %                 util_mask_temp(util_mask_temp ~= class) = 0;
    %                 figure; imagesc(util_mask_temp);
    %                 title(['Drone Image - Mask for class ' num2str(class)]);

                    dtwResult2 = dtwResult2 + dtw(pdfs{find(droneClasses == class, 1)}, normY, 100);
                    droneClassesTemp(droneClassesTemp == class) = [];
                else
                    dtwResult2 = dtwResult2 + dtw(zeros(size(normY)), normY, 100);
                end
                dbg = 1;
            end
        end
        
        if ~isempty(droneClassesTemp)
            for cla2 = 1:length(droneClassesTemp)
                dtwResult2 = dtwResult2 + dtw(zeros(size(normY)), pdfs{find(droneClasses == droneClassesTemp(cla2), 1)}, 100);
            end
        end
                
        dtwResult3 = dtw(finalPdfNorm, finalPdfNorm.*geoFinalPdfNorm, 100);
        dtwResult4 = dtw(finalPdfNorm, geoFinalPdfNorm, 100);
%         figure;
%         surf(x22,x11,geoFinalPdf);
%         title('Final PDF');
% 
%         figure;
%         surf(x22,x11,geoFinalPdfNorm);
%         title('Normalized Final PDF');
        
        
%         BW_Target = edge(geo_mask, 'canny', [], sqrt(2));
%         
%         BW_Template_mean=BW_Template-mean(mean(BW_Template));
%         BW_Target_mean=BW_Target-mean(mean(BW_Target));

        resultsMatrix{ponY, ponX}.out = [sum(KLDiv(finalPdfNorm, geoFinalPdfNorm)) dtwResult dtwResult2 dtwResult3 dtwResult4];
%         corrMat = xcorr2(BW_Target_mean,BW_Template_mean);
%         figure; mesh(corrMat);
%         resultsMatrix{ponY, ponX}.corr = max(corrMat(:));
        
        save matriz_density_32.mat posMatrix resultsMatrix pX pY ponY ponX
        
        pX = pX + 2*dx;
    end
    pY = pY + 2*dy;
end
fprintf('\Tempo total = %f\n\n', toc);
save matriz_density_22.mat posMatrix resultsMatrix pX pY ponY ponX
dbg = 1;

% matriz de resultados
load matriz_density_22.mat
load R.mat;
results = zeros(size(posMatrix));
for k = 1:size(posMatrix, 1)
    for j = 1:size(posMatrix, 2)
        [lat_srp, lon_srp] = pix2latlon(R_copy, posMatrix{k, j}(1), posMatrix{k, j}(2));
        results(k,j) = m_idist(lon_srp, lat_srp, lon(3), lat(3));
    end
end

% 1 vetor só
allOut = [];
for k = 1:size(resultsMatrix, 1)
    for j = 1:size(resultsMatrix, 2)
        %resultsMatrix{k, j}.out = [resultsMatrix{k, j}.out; resultsMatrix{k, j}.corr];
        allOut = [allOut resultsMatrix{k, j}.out];
    end
end

% matriz de features
dtwM = zeros(size(resultsMatrix));
dtwM2 = zeros(size(resultsMatrix));
dtwM3 = zeros(size(resultsMatrix));
dtwM4 = zeros(size(resultsMatrix));
kulM = zeros(size(resultsMatrix));
results = zeros(size(resultsMatrix));
for k = 1:size(resultsMatrix, 1)
    for j = 1:size(resultsMatrix, 2)
        kulM(k, j) = resultsMatrix{k, j}.out(1);
        dtwM(k, j) = resultsMatrix{k, j}.out(2);
        dtwM2(k, j) = resultsMatrix{k, j}.out(3);
        dtwM3(k, j) = resultsMatrix{k, j}.out(4);
        dtwM4(k, j) = resultsMatrix{k, j}.out(5);
    end
end
figure; imagesc(kulM); colorbar;
figure; imagesc(dtwM); colorbar;
figure; imagesc(dtwM2); colorbar;
figure; imagesc(dtwM3); colorbar;
figure; imagesc(dtwM4); colorbar;

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
