function [yoffSet, xoffSet, corrMat, centro, BW_Template, BW_Target]=edges_and_correlation_BJ(image1,image2,n_case,printFigsCorr)
if printFigsCorr
    figure;
    imshowpair(image2,image1,'montage');
end

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
sigmaValue = sqrt(2);
kValue = 0.7;
entropia_Target = entropy(Target);
entropia_Template = entropy(Template);
sigmaValue_Target = sqrt(entropia_Target/2);
sigmaValue_Template = sqrt(entropia_Template/2);

switch n_case
    case 1 % filtro std, sigma sqrt(entropia_Template/2), canny 
        fprintf('Std Filter + Canny \n');
        Template_filt = normalize_img(stdfilt(Template,ones(15)));
        Target_filt = normalize_img(stdfilt(Target,ones(15)));

        BW_Template = edge(Template_filt, filtro, [], sigmaValue_Template);
        BW_Target = edge(Target_filt, filtro, [], sigmaValue_Target);
        
    case 2        
        fprintf('ILS Filter + Canny + kValue\n');
        lambda_Target = sigmaValue_Target/entropia_Target;
        lambda_Template = sigmaValue_Template/entropia_Template;
        iter = 2;
        p = 0.3;
        eps = 0.0001;  

        Template = im2double(Template);
        Template_filt = ILS_LNorm(Template, lambda_Template, p, eps, iter);

        Target = im2double(Target);
        Target_filt = ILS_LNorm(Target, lambda_Target, p, eps, iter);
        
        Template_filt = im2uint8(Template_filt);
        Target_filt = im2uint8(Target_filt);
        
        [thres_Template,sigm_Template] = generate_threshold(Template_filt);
        [thres_Target,sigm_Target] = generate_threshold(Target_filt);
        
        BW_Target = edge(Target_filt, filtro, thres_Target,sigm_Target);
        BW_Template = edge(Template_filt, filtro, thres_Template,sigm_Template);

    case 3 
        fprintf('Canny + kValue + Sigma\n');
        [thres_Template,sigm_Template] = generate_threshold(Template);
        [thres_Target,sigm_Target] = generate_threshold(Target);
        
        BW_Target = edge(Target, filtro, thres_Target,sigm_Target);
        BW_Template = edge(Template, filtro, thres_Template,sigm_Template);
        
    case 4% filtro median, canny
        fprintf('Median Filter + Canny\n');
        
        BW_Template = edge(Template, filtro, [], sigmaValue_Template);
        BW_Target = edge(Target, filtro, [], sigmaValue_Target);
        
    case 5 % filtro median, canny
        fprintf('Median Filter + Canny + kValue\n');
        Template_filt = medfilt2(Template,[5 5]);
        Target_filt = medfilt2(Target,[5 5]);
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

        BW_TemplateAll = [{}];
        BW_TargetAll = [{}];
        for j = 1:length(kValue)
            [~, threshOut1] = edge(Template_filt, filtro, [], sigmaValue);
            [~, threshOut2] = edge(Target_filt, filtro, [], sigmaValue);
            BW_Template = edge(Template_filt, filtro, kValue(j)*threshOut1, sigmaValue);
            BW_Target = edge(Target_filt, filtro, kValue(j)*threshOut2, sigmaValue);
        end
        
    case 8 % filtro entropy, sigma sqrt(2), canny 
        fprintf('Entropy Filter + Canny + kValue\n');
        Template_filt = normalize_img(entropyfilt(Template,ones(9)));
        Target_filt = normalize_img(entropyfilt(Target,ones(9)));

        BW_Template = edge(Template_filt, filtro, [], sigmaValue);
        BW_Target = edge(Target_filt, filtro, [], sigmaValue);
        
    case 9% filtro std
        fprintf('Std Filter\n');
        Template_filt = normalize_img(stdfilt(Template,ones(9)));
        Target_filt = normalize_img(stdfilt(Target,ones(9)));

        BW_Template = Template_filt;
        BW_Target = Target_filt;
        
    case 10 % filtro entropy
        fprintf('Entropy Filter\n');
        Template_filt = normalize_img(entropyfilt(Template,ones(9)));
        Target_filt = normalize_img(entropyfilt(Target,ones(9)));

        BW_Template = Template_filt;
        BW_Target = Target_filt;
end

function [threshold,sigm] = generate_threshold(I)
            C = 0.2;
            threshold_low = max(1,C*(mean2(I)-std2(I)))/255;
            threshold_high = min(254,C*(mean2(I)+std2(I)))/255;
            threshold = [threshold_low threshold_high];
             sigm = mean(threshold)/C;
        %     sigm = threshold_high/C;
end
    
[r1,c1]=size(BW_Target);
[r2,c2]=size(BW_Template);

BW_Template_mean=BW_Template-mean(mean(BW_Template));
BW_Target_mean=BW_Target-mean(mean(BW_Target));

if printFigsCorr
    figure;
    imshowpair(BW_Target_mean,BW_Template_mean,'montage');
end
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

ptCenterX = floor(unique(i) - (r2/2));
ptCenterY = floor(unique(j) - (c2/2));

% hFig = figure(55);
% hAx  = axes;
% imshow(Target,'Parent', hAx);
% imrect(hAx, [xoffSet, yoffSet, size(Template,2), size(Template,1)]);


if all(corrMat <= 0)
    centro = [0 0];
else
    centro = [ptCenterY(1) ptCenterX(1)];
end
end
