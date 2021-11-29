% switch myFeatureExtractor
if printResults
fprintf('\nExtracting features...\n');
end

% Extracting superpixels centroids and adjacency matrix
% 
numRows = size(L,1);
numCols = size(L,2);
superPixels = cell(1, N);
totalAdj = 0;
for labelVal = 1:N
    [y, x] = find(L == labelVal);
    superPixels{labelVal}.centX = round(mean(x));
    superPixels{labelVal}.centY = round(mean(y));
    superPixels{labelVal}.sizeX = max(x)-min(x);
    superPixels{labelVal}.sizeY = max(y)-min(y);
    nX = min(x+1, numCols);
    nY = min(y+1, numRows);
    pX = max(x-1, 1);
    pY = max(y-1, 1);
    adjacents = zeros(length(x), 4);
    adjacents(:, 1) = L(sub2ind(size(L), y, nX));
    adjacents(:, 2) = L(sub2ind(size(L), nY, x));
    adjacents(:, 3) = L(sub2ind(size(L), y, pX));
    adjacents(:, 4) = L(sub2ind(size(L), pY, x));
    uniAdjs = unique(adjacents);
    uniAdjs = uniAdjs(uniAdjs > labelVal);
    totalAdj = totalAdj + length(uniAdjs);
    superPixels{labelVal}.adj = uniAdjs;
end
clear nX nY pX pY x y
% toc
dgb = 1;

pixels = [];

if printResults
tic;
end
% tic;
% Média dos centróides
    if (myFeatureExtractor == 1 || myFeatureExtractor == 4)
%         tic;
        nColors = 3;
        nFts = 5;
        nSpaces = 7;
        pixels = zeros(N, nSpaces*nFts*nColors);
        for t = 1:nSpaces
            switch t
                case 1
                    A = rgbImage;
                case 2
                    A = rgb2hsv(rgbImage);
                case 3
                    A = rgb2lab(rgbImage);
                case 4
                    A = lin2rgb(rgbImage);
                case 5
                    A = rgb2ycbcr(rgbImage);
                case 6
                    A = rgb2ntsc(rgbImage);
                case 7
                    A = rgb2cmyk(rgbImage, 'ori');
%                     cform = makecform('srgb2cmyk');
%                     A = applycform(lin2rgb(rgbImage), cform); 
                    A = A(:,:,1:3);
                case 8
                    A = rgb2cmyk(rgbImage, 'mod'); 
%                     cform = makecform('srgb2cmyk');
%                     A = applycform(lin2rgb(rgbImage), cform); 
                    A = A(:,:,1:3);
            end
%             A(:,:,4:6) = rgb2hsv(rgbImage);
    %         labImage = rgb2lab(rgbImage);
    %         A(:,:,7) = labImage(:,:,1)/100;
    %         A(:,:,8:9) = (labImage(:,:,2:3)+100)/200;
%             A(:,:,7:9) = rgb2lab(rgbImage);
%             A(:,:,10:12) = lin2rgb(rgbImage);
    %         A(:,:,13:15) = rgb2ycbcr(rgbImage);
    %         ntscImage = rgb2ntsc(rgbImage);
    %         A(:,:,16) = ntscImage(:,:,1);
    %         A(:,:,17) = (ntscImage(:,:,2)+0.5959)/(2*0.5959);
    %         A(:,:,18) = (ntscImage(:,:,3)+0.5229)/(2*0.5229);
    %         A(:,:,16:18) = rgb2ntsc(rgbImage);
%             cform = makecform('srgb2cmyk');
%             cymA = applycform(rgbImage, cform); 
%             A(:,:,19:21) = cymA(:,:,1:3);
            clear cymA labImage ntscImage cform
    %         A(:,:,22) = rgb2gray(rgbImage); % is present in YIQ (NTSC)
%             pixels = zeros(N, size(A, 3)*4);
            numRows = size(A,1);
            numCols = size(A,2);
            gradA = gradient(A);
            
            for labelVal = 1:N
                pTemp = [];
                Idx = idx{labelVal} + (0:size(A, 3)-1)*numRows*numCols;
                pTemp = [pTemp mean(A(Idx))];
                pTemp = [pTemp std(A(Idx))];
                pTemp = [pTemp median(A(Idx))];
                pTemp = [pTemp mean(A(Idx).^2)];
                pTemp = [pTemp mean(gradA(Idx))];
%                 pixels(labelVal, :) = pTemp;
                pixels(labelVal, ((t-1)*(nFts*nColors))+(1:(nFts*nColors))) = pTemp;
            end
            clear A gradA pTemp Idx
        end
%         toc;
    end
    
pixelsOwn = pixels;

% LBP
    if (myFeatureExtractor == 2)
%         tic;
        % maior superpixel
        maiorCol = 0;
        maiorRow = 0;
        spValues = [];
        for labelVal = 1:N
            [row,col] = find(L == labelVal);
            if max(row)-min(row) > maiorRow
                maiorRow = max(row)-min(row);
            end
            if max(col)-min(col) > maiorCol
                maiorCol = max(col)-min(col);
            end
            spValues = [spValues; max(row), min(row), max(row)-min(row),...
                        max(col), min(col), max(col)-min(col)];
        end

        % segmenta e aplica
%         A = rgb2gray(rgbImage);
        A = rgbImage;
        A2 = rgb2gray(rgbImage);
        numRows = size(A,1);
        numCols = size(A,2);
        pixels = [];
        offset20x = (round(0.2*numRows));
        offset20y = (round(0.2*numCols));
        imageL = zeros(numRows+2*offset20x, numCols+2*offset20y);
        
        finalIm = zeros(maiorRow, maiorCol);
        temp = zeros(size(A, 1), size(A, 2), 4);
        tempR = A(:,:,1);
        tempG = A(:,:,2);
        tempB = A(:,:,3);
        tempGray = A2;
%             tempR(L ~= labelVal) = 0;
%             tempG(L ~= labelVal) = 0;
%             tempB(L ~= labelVal) = 0;
%             tempGray(L ~= labelVal) = 0;
        temp(:,:,1) = tempR;
        temp(:,:,2) = tempG;
        temp(:,:,3) = tempB;
        temp(:,:,4) = tempGray;
        
        for labelVal = 1:N
            % aqui eu replico a imagem
            [row,col] = find(L == labelVal);
            currImg = temp(min(row):max(row), min(col):max(col), :);
%             finalIm = padarray(currImg,[max(round((maiorRow/2)-(size(currImg, 1)/2)), 0) max(round((maiorCol/2)-(size(currImg, 2)/2)), 0)], 'circular');
            finalIm = currImg;
            % Aqui eu deixo a imagem centralizada
%             imageL(offset20x:numRows+offset20x-1, offset20y:numCols+offset20y-1) = temp;
%             % figure; imshow(imageL);
%             offsetImx = maiorRow - spValues(labelVal, 3);
%             offsetImy = maiorCol - spValues(labelVal, 6);
%             if mod(offsetImx/2, 1) > 0
%                 offx1 = offsetImx/2-0.5; offx2 = offsetImx/2+0.5;
%             else
%                 offx1 = offsetImx/2; offx2 = offsetImx/2;
%             end
%             if mod(offsetImy/2, 1) > 0
%                 offy1 = offsetImy/2-0.5; offy2 = offsetImy/2+0.5;
%             else
%                 offy1 = offsetImy/2; offy2 = offsetImy/2;
%             end
%             distX = offset20x+spValues(labelVal, 2)-offx1;
%             distY = offset20y+spValues(labelVal, 5)-offy1;
%             finalIm = imageL(distX:distX+maiorRow, distY:distY+maiorCol);
            % figure; imshow(finalIm);
            corTemp = [];
            for cor = 1:4
                corTemp = [corTemp extractLBPFeatures(finalIm(:,:,cor))];
            end
            pixels = [pixels; corTemp];
%             dgb = 1;
        end
        toc;
        pixels = double(pixels);
        dgb = 1;
    end
    
    if (myFeatureExtractor == 4)
        pixelsTemp = pixels;
    end

% LM Filters
    if (myFeatureExtractor == 3 || myFeatureExtractor == 4)
        
        pixels = zeros(N, 5*3*18); % 5 caracteristicas 3 canais 18 filtros
        pixelsAdj = zeros(totalAdj, 2+5*3*18); % 5 caracteristicas 3 canais 18 filtros
        A = rgbImage;
        numRows = size(A,1);
        numCols = size(A,2);
%             tic;
            F = makeLMfilters;
%             myAngles = unique(angF);
%             pixels = zeros(N, 3*4*size(F, 3));
            angCnt = 0;
            colorCnt = 0;
            for col = 1:3
                filtCnt = 0;
                for fil = 1:size(F, 3)
                    fA = imfilter(A(:, :, col), F(:,:,fil));
                    if (angCnt == 0), fAold = fA; end

                    if (angCnt < 5 && fil < 37)
                        fAold = max(fA, fAold);
                        angCnt = angCnt+1;
                        continue;
                    elseif (angCnt == 5 && fil < 37)
                        fA = max(fA, fAold);
                        angCnt = 0;
                    end

                    gradfA = gradient(fA);
%                     disp(1+colorCnt+filtCnt);
                    adjCnt = 1;
                    for labelVal = 1:N
                        colorIdx = idx{labelVal};
                        pixels(labelVal, 1+colorCnt+filtCnt) = mean(fA(colorIdx));
                        pixels(labelVal, 2+colorCnt+filtCnt) = std(fA(colorIdx));
                        pixels(labelVal, 3+colorCnt+filtCnt) = median(fA(colorIdx));
                        pixels(labelVal, 4+colorCnt+filtCnt) = mean(fA(colorIdx).^2);
                        pixels(labelVal, 5+colorCnt+filtCnt) = mean(gradfA(colorIdx));
                        for adj = 1:length(superPixels{labelVal}.adj)
                            currAdj = superPixels{labelVal}.adj(adj);
                            xi = min(superPixels{labelVal}.centX, superPixels{currAdj}.centX);
                            xf = max(superPixels{labelVal}.centX, superPixels{currAdj}.centX);
                            yi = min(superPixels{labelVal}.centY, superPixels{currAdj}.centY);
                            yf = max(superPixels{labelVal}.centY, superPixels{currAdj}.centY);
                            faAdj = fA(yi:yf, xi:xf);
                            gradfAAdj = gradient(faAdj);
                            pixelsAdj(adjCnt, 1) = labelVal;
                            pixelsAdj(adjCnt, 2) = currAdj;
                            pixelsAdj(adjCnt, 3+colorCnt+filtCnt) = mean(faAdj(:));
                            pixelsAdj(adjCnt, 4+colorCnt+filtCnt) = std(faAdj(:));
                            pixelsAdj(adjCnt, 5+colorCnt+filtCnt) = median(faAdj(:));
                            pixelsAdj(adjCnt, 6+colorCnt+filtCnt) = mean(faAdj(:).^2);
                            pixelsAdj(adjCnt, 7+colorCnt+filtCnt) = mean(gradfAAdj(:));
                            adjCnt = adjCnt + 1;
                            dbg = 1;
                        end
                    end
                    
                    filtCnt = filtCnt+5;
                end
                colorCnt = colorCnt+filtCnt;
            end
            
            clear A fA fAold gradfA F faAdj gradfAAdj colorIdx
                
%                 for labelVal = 1:N
%                     [row,col] = find(L == labelVal);
%                     sizeX = max(row)-min(row);
%                     sizeY = max(col)-min(col);
%                     temp = rgbImage;
%                     temp(L ~= labelVal) = 0;
%                     for j=1:3
%                         currImg = temp(min(row):max(row), min(col):max(col), j);
%                         responses = zeros(size(currImg, 1)+size(F, 3), size(currImg, 2)+size(F, 3), size(F, 3));
%                         stdValues = zeros(1, size(F, 3));
%                         entValues = zeros(1, size(F, 3));
%                         for i=1:size(F, 3)
%                             responses(:,:,i) = conv2(currImg, F(:,:,i));
%                             stdValues(i) = std2(responses(:,:,i));
%                             entValues(i) = entropy(responses(:,:,i));
%                         end
%                         meanValues = zeros(1, size(F, 3));
%                         meanValues(1:end) = mean(mean(responses));
%                         medianValues = zeros(1, size(F, 3));
%                         medianValues(1:end) = median(median(responses));
%                         attributes = [meanValues; stdValues; medianValues; entValues];
%                         pixels(labelVal, (j-1)*4*size(F, 3)+1:j*4*size(F, 3)) = reshape(attributes, 1, 4*size(F, 3));
%                     end
%                 end
%             end
%             toc;
%         end
        dgb = 1;
    end
    
    if (myFeatureExtractor == 4)
        pixels = [pixelsTemp pixels];
    end
    
    % AUTOENCODER
    if (myFeatureExtractor == 6)
        if m == 1
            load feats.mat
        else
            load(['feat2_m' num2str(m) '.mat']);
        end
       pixels = feat2;
       clear feat1 feat2
    end
    if (myFeatureExtractor == 5)
%         tic;
        % maior superpixel
        maiorCol = 0;
        maiorRow = 0;
        spValues = [];
        for labelVal = 1:N
            [row,col] = find(L == labelVal);
            if max(row)-min(row) > maiorRow
                maiorRow = max(row)-min(row);
            end
            if max(col)-min(col) > maiorCol
                maiorCol = max(col)-min(col);
            end
            spValues = [spValues; max(row), min(row), max(row)-min(row),...
                        max(col), min(col), max(col)-min(col)];
        end

        % segmenta e aplica
%         A = rgb2gray(rgbImage);
        A = rgbImage;
        A2 = rgb2gray(rgbImage);
        numRows = size(A,1);
        numCols = size(A,2);
        pixels = [];
        offset20x = (round(0.2*numRows));
        offset20y = (round(0.2*numCols));
        imageL = zeros(numRows+2*offset20x, numCols+2*offset20y);
        allSuperpixels = [{}];
        for labelVal = 1:N
            finalIm = zeros(maiorRow, maiorCol);
            temp = A2;
            temp(L ~= labelVal) = 0;
            % aqui eu replico a imagem
            [row,col] = find(L == labelVal);
            currImg = temp(min(row):max(row), min(col):max(col), :);
            finalIm = padarray(currImg,[max(ceil((maiorRow/2)-(size(currImg, 1)/2)), 0) max(ceil((maiorCol/2)-(size(currImg, 2)/2)), 0)], 'circular');
            finalIm = finalIm(1:maiorRow, 1:maiorCol);
%             finalIm = currImg;
            allSuperpixels = [allSuperpixels; finalIm];

        end
        clear A A2
        hiddenSize1 = 100;
        autoenc1 = trainAutoencoder(allSuperpixels,hiddenSize1, ...
            'MaxEpochs',100, ...
            'L2WeightRegularization',0.004, ...
            'SparsityRegularization',4, ...
            'SparsityProportion',0.15, ...
            'ScaleData', false);
        feat1 = encode(autoenc1,allSuperpixels)';
        save(['feat1_m' num2str(m) '.mat'], 'feat1');
        feat1 = feat1';
        hiddenSize2 = 50;
        autoenc2 = trainAutoencoder(feat1,hiddenSize2, ...
            'MaxEpochs',100, ...
            'L2WeightRegularization',0.002, ...
            'SparsityRegularization',4, ...
            'SparsityProportion',0.1, ...
            'ScaleData', false);
        feat2 = encode(autoenc2,feat1)';
        save(['feat2_m' num2str(m) '.mat'], 'feat2');
%         softnet = trainSoftmaxLayer(feat2,tTrain,'MaxEpochs',400);
%         stackednet = stack(autoenc1,autoenc2,softnet);
        toc;
        pixels = feat1';
        dgb = 1;
    end
    
    if (myFeatureExtractor == 7)

%         tic;
        % maior superpixel
        maiorCol = 0;
        maiorRow = 0;
        spValues = [];
        for labelVal = 1:N
            [row,col] = find(L == labelVal);
            if max(row)-min(row) > maiorRow
                maiorRow = max(row)-min(row);
            end
            if max(col)-min(col) > maiorCol
                maiorCol = max(col)-min(col);
            end
            spValues = [spValues; max(row), min(row), max(row)-min(row),...
                        max(col), min(col), max(col)-min(col)];
        end

        % segmenta e aplica
        A = rgbImage;
        A2 = rgb2gray(rgbImage);
        numRows = size(A,1);
        numCols = size(A,2);
        pixels = [];
        offset20x = (round(0.2*numRows));
        offset20y = (round(0.2*numCols));
        imageL = zeros(numRows+2*offset20x, numCols+2*offset20y);
        
        finalIm = zeros(maiorRow, maiorCol);
        temp = zeros(size(A, 1), size(A, 2), 4);
        tempR = A(:,:,1);
        tempG = A(:,:,2);
        tempB = A(:,:,3);
        tempGray = A2;
%             tempR(L ~= labelVal) = 0;
%             tempG(L ~= labelVal) = 0;
%             tempB(L ~= labelVal) = 0;
%             tempGray(L ~= labelVal) = 0;
        temp(:,:,1) = tempR;
        temp(:,:,2) = tempG;
        temp(:,:,3) = tempB;
        temp(:,:,4) = tempGray;
        
        for labelVal = 1:N            
            % aqui eu replico a imagem
            [row,col] = find(L == labelVal);
            currImg = temp(min(row):max(row), min(col):max(col), :);
%             finalIm = padarray(currImg,[max(round((maiorRow/2)-(size(currImg, 1)/2)), 0) max(round((maiorCol/2)-(size(currImg, 2)/2)), 0)], 'circular');
            finalIm = currImg;
            %         glcms = graycomatrix(rgbImage);
            offsets = [0 1; -1 1;-1 0;-1 -1];
            [glcms,SI] = graycomatrix(finalIm(:,:,4),'Offset',offsets);
            [out] = GLCM_Features4(glcms,0);
            fields = fieldnames(out);
            corTemp = [];
            for i=1:length(fields)
                corTemp = [corTemp getfield(out,fields{i})];
            end
            pixels = [pixels; corTemp];
%             dgb = 1;
        end
        toc;
        pixels = double(pixels);
        pixels(isnan(pixels)) = 0;
        dgb = 1;
    end
    
    if (myFeatureExtractor == 4)
        pixelsTemp = pixels;
    end
    
    if (myNormalization == 1 || myNormalization == 3)
        % por vetor de caracteristicas (por coluna)
        if (m == 1)
            ftMean = mean(pixels);
            ftStd = std(pixels);
            parameters.featMean = ftMean;
            parameters.featStd = ftStd;
        end
        pixels = (pixels-parameters.featMean)./(2.*(3.*parameters.featStd+1)); 
        pixels(isnan(pixels)) = 0;
        clear ftMean ftStd
    end
    
    if (myNormalization == 2 || myNormalization == 3)
        % por superpixel (por linha)
        pixels = (pixels-mean(pixels, 2))./(2.*(3.*std(pixels, 0, 2)+1));
    end
    
    clear pixelsTemp
       
    pixels(isnan(pixels)) = 0;
    if printResults
    fprintf('Execution time for feature extraction: %f s\n', toc);
    end
    dbg = 1;
%         teste = rgbImage(:, :, 1);
% end