%      close all;
 fprintf('\nImage %d:\n', m)
 
 fprintf('\nReading image...\n');
 tic;

 % obtem a imagem RGB, altera seu tamanho e mostra a imagem original (rows x cols x bands)
 currentImage = strtrim(imageNames{m});   
 rgbImage = imread(currentImage); % im2double() - converte pixels para double
%  rgbImage = imresize(rgbImage,0.25); % diminui o tamanho da imagem para diminuir os calculos
 originalRgbImage = rgbImage;
%  [L,Centers] = imsegkmeans(rgbImage,7);
%  B = labeloverlay(rgbImage,L);
%  figure; imshow(B);
%   wavelength = 2.^(0:5) * 3;
%  orientation = 0:45:135;
%  g = gabor(wavelength,orientation);
%  I = rgb2gray(im2single(rgbImage));
%  gabormag = imgaborfilt(I,g);
% %  montage(gabormag,'Size',[4 6]);
%  for i = 1:length(g)
%     sigma = 0.5*g(i).Wavelength;
%     gabormag(:,:,i) = imgaussfilt(gabormag(:,:,i),3*sigma); 
%  end
% %  montage(gabormag,'Size',[4 6])
%  nrows = size(rgbImage,1);
% ncols = size(rgbImage,2);
% [X,Y] = meshgrid(1:ncols,1:nrows);
% featureSet = cat(3,I,gabormag,X,Y);
% L2 = imsegkmeans(featureSet,7,'NormalizeInput',true);
% C = labeloverlay(rgbImage,L2);
% figure; imshow(C);
 rgbImage = im2double(rgbImage);
 if (useCropped && m ~=  1)
     rgbImage2 = rgbImage;
     [~, currImgName] = fileparts(currentImage);
     currentImageCrop = ['./Images/Test/Cropped/ComEqualizacao/cropped_image_' pad(num2str(m), 2, 'left', '0') '_' currImgName '.JPG'];
     rgbImageTempCrop = im2double(imread(currentImageCrop));
     rgbImageCrop = imresize(rgbImageTempCrop, size(rgbImage(:,:,1)));
     rgbImage = rgbImageCrop;
     rgbImageCrop = rgbImage2;
     originalRgbImage = rgbImage;
     clear rgbImage2;
 end
 fprintf('Execution time for reading image: %f s\n', toc);
 % abaixo, plota a imagem original
%  if plotsIM
%      figure('Renderer', 'painters', 'Position', [100 100 800 700]); set(gcf,'color','w');
% %          annotation('textbox', [0.4, 0.88, 0.25, 0.1], 'String', ['IMAGE ' num2str(m)], 'BackgroundColor', [1, 1, 1], 'FontSize', 20, 'HorizontalAlignment', 'center');
% %      subplot(2,4,1);
%      imshow(rgbImage);
%      title('Original Image')
%  end