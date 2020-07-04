isHSV = 0;
isLab = 0;
fprintf('\nChanging color space...\n');
tic;
switch myColorSpace
    case 0
        fprintf('Color space did not change, using RGB\n');
        % do nothing
    case 1
        
         rgbImage = rgb2hsv(rgbImage);
         isHSV = 1;
%         if plotsIM
%             subplot(2,4,2);
%             imshow(rgbImage,'InitialMagnification',67);
%             title('HSV Image')
%         end
    case 2
%         fprintf('\nChanging color space...\n');
        rgbImage = rgb2lab(rgbImage);
        isLab = 1;
%         if plotsIM
%             subplot(2,4,2);
%             imshow(rgbImage,'InitialMagnification',67);
%             title('Lab Image')
%         end
    case 3
%         fprintf('\nChanging color space...\n');
        rgbImage = lin2rgb(rgbImage);
        issRGB = 1;
%         if plotsIM
% %             subplot(2,4,2);
%             figure;
%             imshow(rgbImage,'InitialMagnification',67);
%             title('sRGB Image')
%         end
end
fprintf('Execution time for changing color space: %f s\n', toc);