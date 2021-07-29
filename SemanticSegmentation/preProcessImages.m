if (myPreprocess == 1 || myPreprocess == 2)
    if printResults
    fprintf('\nPreprocessing image...\n');
    tic;
    end
    % Histogram equalization
    if (myPreprocess == 2)
        labImage = rgb2lab(rgbImage);
        % [counts,binLocations] = imhist(rgbImage);
        % figure; bar(binLocations,counts);
        L = labImage(:,:,1)/100;
        L = adapthisteq(L);
        labImage(:,:,1) = L*100;
        rgbImage = lab2rgb(labImage);
        % figure; montage({rgbImage, rgbImage2});
        % [counts,binLocations] = imhist(rgbImage2);
        % figure; bar(binLocations,counts);
        clear labImage L
        
%         if exist('geoImg','var')
%             labImage = rgb2lab(geoImg);
%             % [counts,binLocations] = imhist(rgbImage);
%             % figure; bar(binLocations,counts);
%             L = labImage(:,:,1)/100;
%             L = adapthisteq(L);
%             labImage(:,:,1) = L*100;
%             geoImg = lab2rgb(labImage);
%             % figure; montage({rgbImage, rgbImage2});
%             % [counts,binLocations] = imhist(rgbImage2);
%             % figure; bar(binLocations,counts);
%             clear labImage
%             oriEqImg = geoImg;
%             clear geoImg
%         end
    end

    if (m == 1)
        % Compute histogram of the reference image
        hgram = zeros(3, 256);
        for currentChan = 1:3
            hgram(currentChan, :) = imhist(rgbImage(:,:,currentChan), 256);
%             [~, hgram(currentChan, :)] = imhist(oriEqImg(:,:,currentChan), 256);
        end
        parameters.eqHist = hgram;
%         oriEqImg = rgbImage;
%         parameters.eqImg = oriEqImg;
%         clear oriEqImg
    else
        % Histogram matching (TODO: use cumuative density function)
    %     [hm, im3] = homography(oriEqImg, rgbImage);
    
    
    % parte para comentar cropped
        for k = 1:size(rgbImage, 3) % Process one color channel at a time
            hgramToUse = k;
            % Use A to store output, to save memory
            rgbImage(:, :, k) = histeq(rgbImage(:,:,k), parameters.eqHist(hgramToUse,:));
        end
    % fim da parte para comentar cropped
        
        
%         rgbImage = imhistmatch(rgbImage, oriEqImg, 256);
    %     figure; montage({rgbImage, rgbImage3});
    %     [counts,binLocations] = imhist(rgbImage3);
    %     figure; bar(binLocations,counts);
    end
    if printResults
    fprintf('Execution time for preprocessing image: %f s\n', toc);
    end
end
