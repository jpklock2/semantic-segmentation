switch myFilter
    % Original
    case 1
        % do nothing
        
    case 3 % Bilateral
        labImage = rgb2lab(rgbImage);
        minRow = std(labImage(:,:,1), [], 2);
        minRowIdx = find(minRow == min(minRow), 1);
        minCol = std(labImage(:,:,1));
        minColIdx = find(minCol == min(minCol), 1);
        square = [max(1, minColIdx-25) max(1, minRowIdx-25) min(size(rgbImage, 2), minColIdx+25) min(size(rgbImage, 1), minRowIdx+25)];
        squareCrop = [square(1) square(2) square(3)-square(1) square(4)-square(2)];
        patch = imcrop(labImage,squareCrop);
        patchSq = patch.^2;
        edist = sqrt(sum(patchSq,3));
        patchVar = std2(edist).^2;
        DoS2 = 4*patchVar;
        smoothedLAB = imbilatfilt(labImage,DoS2,4);
        smoothedRBG = lab2rgb(smoothedLAB,'Out','double');
        clear smoothedLAB smoothedRBG labImage
 
    case 2 % Kuwahara
        hsvImage = rgb2hsv(rgbImage);
        k = 2; % 5, 9, 13, ... = (4*k+1)
        kuwaImg = Kuwahara(hsvImage(:,:,3),4*k+1);
        hsvImage(:,:,3) = kuwaImg;
        kuwaImg = hsv2rgb(hsvImage);
        rgbImage = kuwaImg;
        clear hsvImage kuwaImg
 
    case 5 % Anisiotropica
         hsvImage = rgb2hsv(rgbImage);
         diffImg = imdiffusefilt(hsvImage(:,:,3),'NumberOfIterations',15);
         hsvImage(:,:,3) = diffImg;
         diffImg = hsv2rgb(hsvImage);
         rgbImage = diffImg;
         clear diffImg hsvImage
        
    case 4 %ILS_Welsch
        entropia_img = entropy(rgbImage);
        minRow = std(im2double(rgbImage(:,:,1)), [], 2);
        minCol = std(im2double(rgbImage(:,:,1)));
        
        lambda = entropy(rgbImage)/4;
        gamma = mean([max(minRow) max(minCol)]);
        iter = 2;
        
%         lsqImage = ILS_LNorm(im2double(rgbImage), lambda, p, eps, iter);
        lsqImage = ILS_Welsch(im2double(rgbImage), lambda, gamma, iter);

        rgbImage = im2uint8(lsqImage);
        clear lsqImage
end
fprintf('Execution time for applying filter: %f s\n', toc);
