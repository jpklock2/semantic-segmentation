
nP = 10;
posMatrix = cell(nP);
pY = cropSize(1)+dy;
for ponY = 1:nP
    pX = cropSize(2)+dx;
    for ponX = 1:nP
        pX = pX + 2*dx;
        posMatrix{ponY, ponX} = [pY pX];
    end
    pY = pY + 2*dy;
end



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
%         [newGeoAdjacencies, ~, ~, maskGeoCrop, classesGeoCrop] = getAdjacencies(subImg, parameters, filter);
        title(['Figure X = ' num2str(ponX) ', Y = ' num2str(ponY)]);
%         spCent = zeros(length(newGeoAdjacencies), 2);
%         sizesSpGeo = zeros(1, length(newGeoAdjacencies));
%         oldGeoAdjacencies = newGeoAdjacencies;
%         for adj = 1:length(newGeoAdjacencies)
%             newGeoAdjacencies{adj}.centX = newGeoAdjacencies{adj}.centX + x0-1;
%             newGeoAdjacencies{adj}.centY = newGeoAdjacencies{adj}.centY + y0-1;
%             spCent(adj, :) = [newGeoAdjacencies{adj}.centY newGeoAdjacencies{adj}.centX];
%             sizesSpGeo(adj) = newGeoAdjacencies{adj}.sizeX;
%             newGeoAdjacencies{adj}.cla = classesGeoCrop(adj);
%         end
%         geoAdjacencies = newGeoAdjacencies;
        
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
%         LmaskGeoCrop = zeros(size(maskGeoCrop));
%         for mx = 1:length(classesGeoCrop)
%             LmaskGeoCrop(maskGeoCrop == mx) = classesGeoCrop(mx);
%         end
        
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
%         subImage.rec = recGeoCrop;
%         subImage.mask = maskGeoCrop;
%         subImage.class = classesGeoCrop;
%         subImage.Lmask = LmaskGeoCrop;
        
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
            resultsMatrix{ponY, ponX}.out = [simbin(BW_Template,BW_Target); max(corrMat(:))];
        catch
            dbg = 1;
        end
        
%         figure; mesh(corrMat);
%         resultsMatrix{ponY, ponX}.corr = max(corrMat(:));
        
        save(['matriz_corre_img' num2str(imgCnt) '.mat'], 'posMatrix', 'resultsMatrix', 'pX', 'pY', 'ponY', 'ponX', 'BW_Template');
        
        pX = pX + 2*dx;
    end
    pY = pY + 2*dy;
end
fprintf('\Tempo total = %f\n\n', toc);
dbg = 1;
