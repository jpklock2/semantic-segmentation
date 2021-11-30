function [centro1, centro2, centro3] = edges_and_correlation_v3(image1, image2, image3, imageMap)

Target = imageMap;
Template1 = image1;
Template2 = image2;
Template3 = image3;

function [img_out] = normalize_img(img_in)
    max_img = max(img_in(:));
    min_img = min(img_in(:));
    img_out = (img_in - min_img) / (max_img - min_img);
end

% filtro = 'canny';
% kValue = 0.7;
% sigmaValue = sqrt(2);
% entropia_Target = entropy(Target);
% entropia_Template = entropy(Template1);
% sigmaValue_Target = sqrt(entropia_Target/2);
% sigmaValue_Template = sqrt(entropia_Template/2);


BW_Template1 = edge(Template1, 'canny', [], sqrt(2));
BW_Template2 = edge(Template2, 'canny', [], sqrt(2));
BW_Template3 = edge(Template3, 'canny', [], sqrt(2));
BW_Target = edge(Target, 'canny', [], sqrt(2));

[r1,c1] = size(BW_Target);
[r21,c21] = size(BW_Template1);
[r22,c22] = size(BW_Template2);
[r23,c23] = size(BW_Template3);

BW_Template_mean1 = BW_Template1 - mean(mean(BW_Template1));
BW_Template_mean2 = BW_Template2 - mean(mean(BW_Template2));
BW_Template_mean3 = BW_Template3 - mean(mean(BW_Template3));
BW_Target_mean = BW_Target - mean(mean(BW_Target));

% figure (81);
% imshowpair(BW_Target_mean,BW_Template_mean,'montage');
    
% corrMat = normxcorr2(image22,Nimage1);

% if exist('redArea', 'var')
%     [sizeY, sizeX] = size(BW_Target_mean);
%     mapY = redLimits(2):redLimits(4);
%     mapX = redLimits(1):redLimits(3);
%     cropY = cropSize(2):cropSize(4);
%     cropX = cropSize(1):cropSize(3);
%     offYi = mapY(1)-cropY(1);
%     offYf = cropY(end)-mapY(end);
%     offXi = mapX(1)-cropX(1);
%     offXf = cropX(end)-mapX(end);
%     BW_Target_mean = BW_Target_mean(offYi+1:sizeY-offYf, offXi+1:sizeX-offXf);
% end

corrMat1 = xcorr2(BW_Target_mean, BW_Template_mean1);
corrMat2 = xcorr2(BW_Target_mean, BW_Template_mean2);
corrMat3 = xcorr2(BW_Target_mean, BW_Template_mean3);

cropCorrMat1 = corrMat1(r21/2:end-r21/2, c21/2:end-c21/2);
cropCorrMat2 = corrMat2(r22/2:end-r22/2, c22/2:end-c22/2);
cropCorrMat3 = corrMat3(r23/2:end-r23/2, c23/2:end-c23/2);

% figure(9)
% mesh(corrMat);

% [r3, c3] = find(corrMat==max(corrMat(:)));
% ypeak=r3;
% xpeak=c3;

% yoffSet = ypeak-size(Template1,1);
% xoffSet = xpeak-size(Template1,2);

% [r,c]=max(corrMat);

% i=c(c3);
% j=c3;

% [ssr,snd] = max(corrMat(:));
% [ij,ji] = ind2sub(size(corrMat),snd);

% ptCenterX = floor(i - (r21/2));
% ptCenterY = floor(j - (c21/2));

% centro_old = [ptCenterX, ptCenterY];

% if exist('redArea', 'var')
%     redCorrMat = cropCorrMat.*redArea;
%     [ptCenterX, ptCenterY] = find(redCorrMat == max(redCorrMat(:)));
% end

[ptCenterX1, ptCenterY1] = find(cropCorrMat1 == max(cropCorrMat1(:)));
[ptCenterX2, ptCenterY2] = find(cropCorrMat2 == max(cropCorrMat2(:)));
[ptCenterX3, ptCenterY3] = find(cropCorrMat3 == max(cropCorrMat3(:)));

% hFig = figure(55);
% hAx  = axes;
% imshow(Target,'Parent', hAx);
% imrect(hAx, [xoffSet, yoffSet, size(Template,2), size(Template,1)]);

centro1 = [ptCenterX1 ptCenterY1];
centro2 = [ptCenterX2 ptCenterY2];
centro3 = [ptCenterX3 ptCenterY3];

end
