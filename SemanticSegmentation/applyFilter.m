if printResults
fprintf('\nApplying filter...\n');
tic;
end
switch myFilter

     % Original
     case 1
         % do nothing
         
 % Kuwahara
     case 2
         hsvImage = rgb2hsv(rgbImage);
%              grayImage = rgb2gray(rgbImage);
        k = 2; % 5, 9, 13, ... = (4*k+1)
         kuwaImg = Kuwahara(hsvImage(:,:,3),4*k+1);
         hsvImage(:,:,3) = kuwaImg;
         kuwaImg = hsv2rgb(hsvImage);
%          disp(['The SSIM value of the noisy image is ',num2str(ssim(kuwaImg,rgbImage))]);
%             figure; montage({rgbImage,kuwaImg});
%              figure; imshow(kuwaImg);
        rgbImage = kuwaImg;
        clear hsvImage kuwaImg
        
 % Bilateral
     case 3
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
        if DoS2 == 0
            DoS2 = 1;
        end
        smoothedLAB = imbilatfilt(labImage,DoS2,4);
        smoothedRBG = lab2rgb(smoothedLAB,'Out','double');
%         disp(['The SSIM value of the noisy image is ',num2str(ssim(smoothedRBG,rgbImage))]);
%             figure; montage({rgbImage,smoothedRBG});
        rgbImage = smoothedRBG;
        clear smoothedLAB smoothedRBG labImage
        
    case 4
    % Mínimos Quadrados
%         rgbTemp = rgbImage;
%         lambda = 3;
%         iter = 8;
%         p = 0.5;
%         eps = 0.0001;
%         lsqImage = double(ILS_LNorm(rgbImage, lambda, p, eps, iter));
% %         figure; montage({rgbImage, lsqImage, kuwaImg, diffImg});
%         rgbImage = lsqImage;
%         clear lsqImage
%         figure; montage({rgbTemp,rgbImage});
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

 % Anisiotropica
     case 5
%              grayImage = rgb2gray(rgbImage);
        hsvImage = rgb2hsv(rgbImage);
%         [gradThresh,numIter] = imdiffuseest(hsvImage(:,:,3));
%         gradThresh = [0.1098 0.0863 0.0706 0.0627 0.0549];
%         numIter = 5;
%         diffImg = imdiffusefilt(hsvImage(:,:,3), ...
%         'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
        diffImg = imdiffusefilt(hsvImage(:,:,3),'NumberOfIterations',15);
%              diffImg = imdiffusefilt(hsvImage(:,:,3));
         hsvImage(:,:,3) = diffImg;
         diffImg = hsv2rgb(hsvImage);
%          disp(['The SSIM value of the noisy image is ',num2str(ssim(diffImg,rgbImage))]);
%             figure; montage({rgbImage,diffImg});
%              figure; imshow(diffImg);
        rgbImage = diffImg;
        clear diffImg hsvImage
        
end
if printResults
fprintf('Execution time for applying filter: %f s\n', toc);
end