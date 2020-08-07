function [yoffSet, xoffSet, corrMat, centro]=edges_and_correlation(image1,image2,mask,centroids,classes,parameters)

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

function [img_out] = rescale_img(img_in)
    max_img = max(img_in(:));
    min_img = min(img_in(:));
    img_out = (img_in - min_img) / (max_img - min_img);
end

filtro = 'canny';

n = 0;

histValues = sum(histc(mask, unique(mask)), 2);
predominance = histValues(2:end)./sum(histValues(2:end));
sigmaValue = sqrt(2) * (((max(predominance)-0.25)/0.75)*2 + 1);
centroidsReconstructed = ((parameters.pcaCoeffs(:,1:parameters.pcaN)*centroids')+parameters.pcaMean')';
textureCentroids = centroidsReconstructed(:, 85:end);
sigC = 1./(1+exp(-mean(textureCentroids, 2)./(3*mean(std(textureCentroids, [], 2)))));
sigC = (sigC+1)/2;
kValue = sum(sigC(unique(classes)).*predominance);

switch n
    case 0
        Template_filt = rescale_img(stdfilt(Template,ones(9)));
        Target_filt = rescale_img(stdfilt(Target,ones(9)));
%         Template = rescale_img(entropyfilt(Template));
%         Target = rescale_img(entropyfilt(Target));
%         figure(5)
%         imshow(Template_filt)
%         figure(6)
%         imshow(Target_filt)
        [~, threshOut1] = edge(Template_filt, filtro, [], sigmaValue);
        [~, threshOut2] = edge(Target_filt, filtro, [], sigmaValue);
        k = 0.70;
        BW_Template = edge(Template_filt, filtro, k*threshOut1, sigmaValue);
        BW_Target = edge(Target_filt, filtro, k*threshOut2, sigmaValue);
        BW_Mask = edge(mask, filtro, k*threshOut1, sigmaValue);
        BW_Template = BW_Template | BW_Mask;
%         figure(7)
%         imshow(BW_Template)
%         figure(8)
%         imshow(BW_Target)

    case 1
        BW_Template = edge(Template, filtro,([]), sqrt(2));
        BW_Target = edge(Target, filtro,([]), sqrt(2));
%         figure(7)
%         imshow(BW_Template)
%         figure(8)
%         imshow(BW_Target)

        
    case 2
        [~, threshOut1] = edge(Template, filtro);
        [~, threshOut2] = edge(Target, filtro);
        BW_Template = edge(Template, filtro,0.75*threshOut1);
        BW_Target = edge(Target, filtro,0.75*threshOut2);
%         figure(7)
%         imshow(BW_Template)
%         figure(8)
%         imshow(BW_Target)
    otherwise
        BW_Template = edge(Template, filtro);
        BW_Target = edge(Target, filtro);
%         figure(7)
%         imshow(BW_Template)
%         figure(8)
%         imshow(BW_Target)
end

[r1,c1]=size(BW_Target);
[r2,c2]=size(BW_Template);

BW_Template_mean=BW_Template-mean(mean(BW_Template));
BW_Target_mean=BW_Target-mean(mean(BW_Target));
% image22=Template;
% Nimage1=Target;

% corrMat = normxcorr2(image22,Nimage1);
corrMat = xcorr2(BW_Target_mean,BW_Template_mean);

% figure(9)
% mesh(corrMat);

%[r3,c3]=max(max(corrMat))
[r3, c3] = find(corrMat==max(corrMat(:)));
ypeak=r3;
xpeak=c3;
yoffSet = ypeak-size(Template,1);
xoffSet = xpeak-size(Template,2);
[r,c]=max(corrMat);

i=c(c3);
j=c3;


[ssr,snd] = max(corrMat(:));
[ij,ji] = ind2sub(size(corrMat),snd);


ptCenterX = floor(i - (r2/2));
ptCenterY = floor(j - (c2/2));

% hFig = figure(55);
% hAx  = axes;
% imshow(Target,'Parent', hAx);
% imrect(hAx, [xoffSet, yoffSet, size(Template,2), size(Template,1)]);

centro = [ptCenterX ptCenterY];
end
