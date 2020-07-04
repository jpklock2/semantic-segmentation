if (myPreprocess == 1 || myPreprocess == 2)
    fprintf('\nPreprocessing image...\n');
    tic;
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
    end

    if (m == 1)
        oriEqImg = rgbImage;
    else
        % Histogram matching (TODO: use cumuative density function)
    %     [hm, im3] = homography(oriEqImg, rgbImage);
        rgbImage = imhistmatch(rgbImage, oriEqImg, 256);
    %     figure; montage({rgbImage, rgbImage3});
    %     [counts,binLocations] = imhist(rgbImage3);
    %     figure; bar(binLocations,counts);
    end
    fprintf('Execution time for preprocessing image: %f s\n', toc);
end