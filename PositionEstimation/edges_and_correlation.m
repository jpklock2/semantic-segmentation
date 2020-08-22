function [yoffSet, xoffSet, corrMat, centro]=edges_and_correlation(image1,image2,n_case,mask,centroids,classes,parameters)

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

% histValues = sum(histc(mask, unique(mask)), 2);
% predominance = histValues(2:end)./sum(histValues(2:end));
% sigmaValue = sqrt(2) * (((max(predominance)-0.25)/0.75)*2 + 1);
% centroidsReconstructed = ((parameters.pcaCoeffs(:,1:parameters.pcaN)*centroids')+parameters.pcaMean')';
% textureCentroids = centroidsReconstructed(:, 85:end);
% sigC = 1./(1+exp(-mean(textureCentroids, 2)./(3*mean(std(textureCentroids, [], 2)))));
% sigC = (sigC+1)/2;
% kValue = sum(sigC(unique(classes)).*predominance);
sigmaValue = sqrt(2);

switch n_case
    case 1 %filtro std, sigma sqrt(2), canny 
        Template_filt = normalize_img(stdfilt(Template,ones(9)));
        Target_filt = normalize_img(stdfilt(Target,ones(9)));
%         Template = normalize_img(entropyfilt(Template));
%         Target = normalize_img(entropyfilt(Target));
%         figure(5)
%         imshow(Template_filt)
%         figure(6)
%         imshow(Target_filt)
        [~, threshOut1] = edge(Template_filt, filtro, [], sigmaValue);
        [~, threshOut2] = edge(Target_filt, filtro, [], sigmaValue);
        k = 0.70;

        BW_Template = edge(Template_filt, filtro, k*threshOut1, sigmaValue);
        BW_Target = edge(Target_filt, filtro, k*threshOut2, sigmaValue);
%         BW_Mask = edge(mask, filtro, k*threshOut1, sigmaValue);
%         BW_Template = BW_Template | BW_Mask;

%         figure(7)
%         imshow(BW_Template)
%         figure(8)
%         imshow(BW_Target)

    case 2 %filtro entropy, sigma sqrt(2), canny 
        Template_filt = normalize_img(entropyfilt(Template,ones(9)));
        Target_filt = normalize_img(entropyfilt(Target,ones(9)));

        [~, threshOut1] = edge(Template_filt, filtro, [], sigmaValue);
        [~, threshOut2] = edge(Target_filt, filtro, [], sigmaValue);
        k = 0.70;

        BW_Template = edge(Template_filt, filtro, k*threshOut1, sigmaValue);
        BW_Target = edge(Target_filt, filtro, k*threshOut2, sigmaValue);
        
    case 3 %filtro std
        Template_filt = normalize_img(stdfilt(Template,ones(9)));
        Target_filt = normalize_img(stdfilt(Target,ones(9)));

        BW_Template = Template_filt;
        BW_Target = Target_filt;
        
    case 4 %filtro entropy
        Template_filt = normalize_img(entropyfilt(Template,ones(9)));
        Target_filt = normalize_img(entropyfilt(Target,ones(9)));

        BW_Template = Template_filt;
        BW_Target = Target_filt;

    case 5 % filtro median, canny
        Template_filt = medfilt2(Template,[3 3]);
        Target_filt = medfilt2(Target,[3 3]);
        [~, threshOut1] = edge(Template_filt, filtro, [], sigmaValue);
        [~, threshOut2] = edge(Target_filt, filtro, [], sigmaValue);
        k = 0.70;

        BW_Template = edge(Template_filt, filtro, k*threshOut1, sigmaValue);
        BW_Target = edge(Target_filt, filtro, k*threshOut2, sigmaValue);

        
    case 6 %filtro median
        Template_filt = medfilt2(Template,[5 5]);
        Target_filt = medfilt2(Target,[5 5]);
        
        BW_Template = Template_filt;
        BW_Target = Target_filt;
        
    otherwise %filtro median, canny
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
corrMat = xcorr2(BW_Target_mean,BW_Template_mean);

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

hFig = figure(55);
hAx  = axes;
imshow(Target,'Parent', hAx);
imrect(hAx, [xoffSet, yoffSet, size(Template,2), size(Template,1)]);

centro = [ptCenterX ptCenterY];
end
