function [yoffSet, xoffSet, corrMat, centro, centro_old] = edges_and_correlation_v2(image1, image2, n_case, redLimits, cropSize, redArea)

if size(image1,3)==3
    image1=rgb2gray(image1);
end
if size(image2,3)==3
    image2=rgb2gray(image2);
end

%Targed = geo image
%Template = uav image
if size(image1)>size(image2)
    Target=image1;
    Template=image2;
else
    Target=image2;
    Template=image1;
end

function [img_out] = normalize_img(img_in)
    max_img = max(img_in(:));
    min_img = min(img_in(:));
    img_out = (img_in - min_img) / (max_img - min_img);
end

filtro = 'canny';
kValue = 0.7;
sigmaValue = sqrt(2);

switch n_case
    case 1 % filtro std, sigma sqrt(2), canny 
        fprintf('Std Filter + Canny + kValue\n');
        Template_filt = normalize_img(stdfilt(Template,ones(9)));
        Target_filt = normalize_img(stdfilt(Target,ones(9)));

        [~, threshOut1] = edge(Template_filt, filtro, [], sigmaValue);
        [~, threshOut2] = edge(Target_filt, filtro, [], sigmaValue);
        
        BW_Template = edge(Template_filt, filtro, kValue*threshOut1, sigmaValue);
        BW_Target = edge(Target_filt, filtro, kValue*threshOut2, sigmaValue);

    case 2 % filtro entropy, sigma sqrt(2), canny 
        fprintf('Entropy Filter + Canny + kValue\n');
        Template_filt = normalize_img(entropyfilt(Template,ones(9)));
        Target_filt = normalize_img(entropyfilt(Target,ones(9)));

        [~, threshOut1] = edge(Template_filt, filtro, [], sigmaValue);
        [~, threshOut2] = edge(Target_filt, filtro, [], sigmaValue);
        
        BW_Template = edge(Template_filt, filtro, kValue*threshOut1, sigmaValue);
        BW_Target = edge(Target_filt, filtro, kValue*threshOut2, sigmaValue);
        
    case 3 % filtro std
        fprintf('Std Filter\n');
        Template_filt = normalize_img(stdfilt(Template,ones(9)));
        Target_filt = normalize_img(stdfilt(Target,ones(9)));

        BW_Template = Template_filt;
        BW_Target = Target_filt;
        
    case 4 % filtro entropy
        fprintf('Entropy Filter\n');
        Template_filt = normalize_img(entropyfilt(Template,ones(9)));
        Target_filt = normalize_img(entropyfilt(Target,ones(9)));

        BW_Template = Template_filt;
        BW_Target = Target_filt;

    case 5 % filtro median, canny
        fprintf('Median Filter + Canny + kValue\n');
        Template_filt = medfilt2(Template,[3 3]);
        Target_filt = medfilt2(Target,[3 3]);
        
        [~, threshOut1] = edge(Template_filt, filtro, [], sigmaValue);
        [~, threshOut2] = edge(Target_filt, filtro, [], sigmaValue);

        BW_Template = edge(Template_filt, filtro, kValue*threshOut1, sigmaValue);
        BW_Target = edge(Target_filt, filtro, kValue*threshOut2, sigmaValue);

        
    case 6 % filtro median
        fprintf('Median Filter\n');
        Template_filt = medfilt2(Template,[5 5]);
        Target_filt = medfilt2(Target,[5 5]);
        
        BW_Template = Template_filt;
        BW_Target = Target_filt;
        
    case 7 % filtro entropy, sigma sqrt(2), canny 
        fprintf('Entropy Filter + Canny + kValue\n');
        Template_filt = normalize_img(entropyfilt(Template,ones(9)));
        Target_filt = normalize_img(entropyfilt(Target,ones(9)));
        
    otherwise % filtro median, canny
        fprintf('Median Filter + Canny\n');
        Template_filt = medfilt2(Template,[5 5]);
        Target_filt = medfilt2(Target,[5 5]);
        
        BW_Template = edge(Template_filt, filtro);
        BW_Target = edge(Target_filt, filtro);

end

[r1,c1]=size(BW_Target);
[r2,c2]=size(BW_Template);

BW_Template_mean=BW_Template-mean(mean(BW_Template));
BW_Target_mean=BW_Target-mean(mean(BW_Target));


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

corrMat = xcorr2(BW_Target_mean, BW_Template_mean);
cropCorrMat = corrMat(r2/2:end-r2/2, c2/2:end-c2/2);

% figure(9)
% mesh(corrMat);

[r3, c3] = find(corrMat==max(corrMat(:)));
ypeak=r3;
xpeak=c3;

yoffSet = ypeak-size(Template,1);
xoffSet = xpeak-size(Template,2);

[r,c]=max(corrMat);

i=c(c3);
j=c3;

% [ssr,snd] = max(corrMat(:));
% [ij,ji] = ind2sub(size(corrMat),snd);

ptCenterX = floor(i - (r2/2));
ptCenterY = floor(j - (c2/2));

centro_old = [ptCenterX, ptCenterY];

if exist('redArea', 'var')
    redCorrMat = cropCorrMat.*redArea;
    [ptCenterX, ptCenterY] = find(redCorrMat == max(redCorrMat(:)));
end



% hFig = figure(55);
% hAx  = axes;
% imshow(Target,'Parent', hAx);
% imrect(hAx, [xoffSet, yoffSet, size(Template,2), size(Template,1)]);

centro = [ptCenterX ptCenterY];
end
