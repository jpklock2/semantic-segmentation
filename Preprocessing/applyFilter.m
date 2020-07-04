fprintf('\nApplying filter...\n');
tic; 
switch myFilter

     % Original
     case 0
         % do nothing

     % Bilateral
     case 1
        labImage = rgb2lab(rgbImage);
        patch = imcrop(labImage,[34,71,60,55]);
        patchSq = patch.^2;
        edist = sqrt(sum(patchSq,3));
        patchVar = std2(edist).^2;
        DoS2 = 4*patchVar;
        smoothedLAB = imbilatfilt(labImage,DoS2,4);
        smoothedRBG = lab2rgb(smoothedLAB,'Out','double');
%         disp(['The SSIM value of the noisy image is ',num2str(ssim(smoothedRBG,rgbImage))]);
%             figure; montage({rgbImage,smoothedRBG});
        rgbImage = smoothedRBG;

 % Kuwahara
     case 2
         hsvImage = rgb2hsv(rgbImage);
%              grayImage = rgb2gray(rgbImage);
        k = 2;
         kuwaImg = Kuwahara(hsvImage(:,:,3),4*k+1);
         hsvImage(:,:,3) = kuwaImg;
         kuwaImg = hsv2rgb(hsvImage);
%          disp(['The SSIM value of the noisy image is ',num2str(ssim(kuwaImg,rgbImage))]);
%             figure; montage({rgbImage,kuwaImg});
%              figure; imshow(kuwaImg);
        rgbImage = kuwaImg;

 % Anisiotropica
     case 3
%              grayImage = rgb2gray(rgbImage);
        hsvImage = rgb2hsv(rgbImage);
%         [gradThresh,numIter] = imdiffuseest(hsvImage(:,:,3),'ConductionMethod','exponential');
        gradThresh = [0.1098 0.0863 0.0706 0.0627 0.0549];
        numIter = 5;
        diffImg = imdiffusefilt(hsvImage(:,:,3),'ConductionMethod','exponential', ...
        'GradientThreshold',gradThresh,'NumberOfIterations',numIter);
%              diffImg = imdiffusefilt(hsvImage(:,:,3));
         hsvImage(:,:,3) = diffImg;
         diffImg = hsv2rgb(hsvImage);
%          disp(['The SSIM value of the noisy image is ',num2str(ssim(diffImg,rgbImage))]);
%             figure; montage({rgbImage,diffImg});
%              figure; imshow(diffImg);
        rgbImage = diffImg;

end
fprintf('Execution time for applying filter: %f s\n', toc);