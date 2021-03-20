function [yoffSet, xoffSet, corrMat, centro, kValue, kValue1, kValue2, kValue3, BW_Template, BW_Target]=edges_and_correlation(image1,image2,n_case,mask,util_mask,rot_mask,ori_mask,ori_mask_geo,centroids,classes,classesGeo,parameters,filter,expe,geo_mask,adjacencies,ftOwn,ftAdj,geoAdjacencies,ftGeoOwn,ftGeoAdj,LmaskGeo,cropSize,utilCropSize,rot_img,geo_img, centroidsGeo)

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

if exist('mask','var')
    
% Sigma baseado em predominancia
histValues = sum(histc(util_mask, unique(util_mask)), 2);
currClasses = unique(util_mask);
if currClasses(1) == 0
    currClasses(1) = [];
    histValues = histValues(2:end);
end
predominance = histValues./sum(histValues);
minProb = 1/length(predominance);
% sigmaValue = sqrt(2) * (((max(predominance)-minProb)/(1-minProb+1e-10))*2 + 1);
sigmaValuePred = ((max(predominance)-minProb)/(1-minProb+1e-10))*2 + 0.5;

% predominancia georreferenciada
histValuesGeo = sum(histc(geo_mask, unique(geo_mask)), 2);
currClassesGeo = unique(geo_mask);
if currClassesGeo(1) == 0
    currClassesGeo(1) = [];
    histValuesGeo = histValuesGeo(2:end);
end
predominanceGeo = histValuesGeo./sum(histValuesGeo);

% Centroides método antigo (média)
% centroidsReconstructed = ((parameters.pcaCoeffs(:,1:parameters.pcaN)*centroids')+parameters.pcaMean')';
% textureCentroids = centroidsReconstructed(:, 85:end);
% sigC = 1./(1+exp(-mean(textureCentroids, 2)./(3*mean(std(textureCentroids, [], 2)))));
% sigC = (sigC+1)/2;
% kValue = sum(sigC(unique(classes)).*predominance);

% Extraindo Gaussianas
inFis = parameters.fis;
cent = zeros(size(inFis.Inputs, 2), size(inFis.Inputs(1).MembershipFunctions, 2));
sig = zeros(size(inFis.Inputs, 2), size(inFis.Inputs(1).MembershipFunctions, 2));
x = -0.5:0.0001:0.5;
mfsInp = {[]};
maxMfsInp = [];
maxMfsInpRed = [];
for j=1:size(inFis.Inputs(1).MembershipFunctions, 2)
    mfs = [];
	for i=1:size(inFis.Inputs, 2)
        mfs = [mfs; gaussmf(x, inFis.Inputs(i).MembershipFunctions(j).Parameters)];
        cent(i,j) = inFis.Inputs(i).MembershipFunctions(j).Parameters(2);
        sig(i,j) = inFis.Inputs(i).MembershipFunctions(j).Parameters(1);
	end
    maxMfsInp = [maxMfsInp; max(mfs)];
    maxMfsInpRed = [maxMfsInpRed; max(mfs(1:2, :))];
    mfsInp = [mfsInp; mfs];
end
mfsInp = mfsInp(2:end);

% Extraindo centroides
% cType = 'bisector';
cType = 'centroid';
allCentroids = [];
allCentroidsRed = [];
for i = 1:size(maxMfsInp, 1)
    allCentroids = [allCentroids; defuzz(x,maxMfsInp(i, :),cType)];
    allCentroidsRed = [allCentroidsRed; defuzz(x,maxMfsInpRed(i, :),cType)];
end

map = 2;

% Mapeamento 1
if map == 1
    sigC = 1./(1+exp(-allCentroids./(std(allCentroids))));
    sigC = (sigC+1)/2;
    sigCRed = 1./(1+exp(-allCentroidsRed./(std(allCentroidsRed))));
    sigCRed = (sigCRed+1)/2;
end

% Mapeamento 2
if map == 2
    sigC = mapminmax(allCentroids', 0.5, 1)';
    sigCRed = mapminmax(allCentroidsRed', 0.5, 1)';
    sigCRed2 = mapminmax(allCentroidsRed', 1, 0.5)';
%     sigKCRed = mapminmax(allCentroidsRed', 0.5, 2)';
    sigKCRed = mapminmax(allCentroidsRed', 2.5, 0.5)';
end

% Método 1: multiplicar as MFs pela predominancia e extrair centroide
% mf = max(predominance.*maxMfsInp(currClasses, :));
% xCentroid1 = defuzz(x,mf,cType);
% idxC = find(x >= xCentroid1);
% kValue1 = mf(idxC(1));
% if map == 1
%     kValue1 = 1./(1+exp(-kValue1./(std(allCentroids))));
%     kValue1 = (kValue1+1)/2;
% elseif map == 2
%     kValue1 = interp1(allCentroids,sigC,kValue1);
% end

% Método 2: multiplicar duas primeiras MFs pela predominancia e extrair centroide
% mf2 = max(predominance.*maxMfsInpRed(currClasses, :));
% xCentroid2 = defuzz(x,mf2,cType);
% idxC2 = find(x >= xCentroid2);
% kValue2 = mf2(idxC2(1));
% if map == 1
%     kValue2 = 1./(1+exp(-kValue2./(std(allCentroidsRed))));
%     kValue2 = (kValue2+1)/2;
% elseif map == 2
%     kValue2 = interp1(allCentroidsRed,sigC,kValue2);
% end

% Atribuindo
% allCentroids = sigC;
allCentroidsRed = sigCRed;
% allCentroidsKRed = sigKCRed;

% % Método 3: multiplicar duas primeiras MFs pela predominancia e extrair centroide
% xCentroid3 = allCentroids(currClasses);
% xCentroid3All = sum(predominance.*xCentroid3);
% % idxC = find(x >= xCentroid3All);
% % kValue3 = mf(idxC(1));
% kValue3 = xCentroid3All;

% Método 4: multiplicar  MFs pela predominancia e extrair centroide
xCentroid4K = allCentroidsRed(currClasses);
xCentroid4KAll = sum(predominance.*xCentroid4K);
% idxC = find(x >= xCentroid4All);
% kValue4 = mf(idxC(1));
kValue4 = xCentroid4KAll;

xCentroid4K1 = sigCRed2(currClasses);
xCentroid4KAll1 = sum(predominance.*xCentroid4K1);
kValue1 = xCentroid4KAll1;

% geo 1
xCentroid4K2 = sigCRed(currClassesGeo);
xCentroid4KAll2 = sum(predominanceGeo.*xCentroid4K2);
kValue2 = xCentroid4KAll2;

% geo 2
xCentroid4K3 = sigCRed2(currClassesGeo);
xCentroid4KAll3 = sum(predominanceGeo.*xCentroid4K3);
kValue3 = xCentroid4KAll3;

% xCentroid4 = allCentroidsKRed(currClasses);
% xCentroid4All = sum(predominance.*xCentroid4);
% sigmaValue4 = xCentroid4All;

% sigmaValue = sigmaValue4;
% sigmaValue = sigmaValuePred;
sigmaValue = sqrt(2);
kValue = kValue4;
% kValue = xCentroid4K(3);
% kValue = 0.7;

minText = find(xCentroid4K == min(xCentroid4K));
minClass = currClasses(minText(1));
minIdx = find(util_mask == minClass, 1);
if isempty(minIdx)
    disp('TA VAZIOOOOOOOO');
end

%% Getting adjacency matrix

% geoAdjacency;
% corrMatrix = geoAdjacency(ftAdj, ftOwn, geoAdjacencies, adjacencies, ftGeoOwn, ftGeoAdj,...
%                           classes, classesGeo, ori_mask, ori_mask_geo, rot_mask, LmaskGeo,...
%                           cropSize, utilCropSize, rot_img, parameters, filter, geo_img);
                      
corrMatrixTM = geoAdjacencyTM_old(ftAdj, ftOwn, geoAdjacencies, adjacencies, ftGeoOwn, ftGeoAdj,...
                          classes, classesGeo, ori_mask, ori_mask_geo, rot_mask, LmaskGeo,...
                          cropSize, utilCropSize, rot_img, parameters, filter, geo_img, util_mask, image1);
                      
% corrMatrixGD = geoAdjacencyGD(ftAdj, ftOwn, geoAdjacencies, adjacencies, ftGeoOwn, ftGeoAdj,...
%                           classes, classesGeo, ori_mask, ori_mask_geo, rot_mask, LmaskGeo,...
%                           cropSize, utilCropSize, rot_img, parameters, filter, geo_img, util_mask, image1, centroidsGeo);

else

sigmaValue = sqrt(2);
kValue = 0.7;

end

switch n_case
    case 1 % filtro std, sigma sqrt(2), canny 
        fprintf('Std Filter + Canny + kValue\n');
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
%         k = 0.70;

        BW_Template = edge(Template_filt, filtro, kValue*threshOut1, sigmaValue);
        BW_Target = edge(Target_filt, filtro, kValue*threshOut2, sigmaValue);
%         BW_Mask = edge(util_mask, filtro, kValue*threshOut1, sigmaValue);
%         BW_Template = BW_Template | BW_Mask;

%         figure(7)
%         imshow(BW_Template)
%         figure(8)
%         imshow(BW_Target)

    case 2 % filtro entropy, sigma sqrt(2), canny 
        fprintf('Entropy Filter + Canny + kValue\n');
        Template_filt = normalize_img(entropyfilt(Template,ones(9)));
        Target_filt = normalize_img(entropyfilt(Target,ones(9)));

%         [~, threshOut1] = edge(Template_filt, filtro, [], sigmaValue);
%         [~, threshOut2] = edge(Target_filt, filtro, [], sigmaValue);
%         k = 0.70;

%         BW_Template = Canny_v0(Template_filt, filter, kValue, expe, minIdx);
%         BW_Target = Canny_v0(Target_filt, filter, kValue, expe, minIdx);

        % CORRETO!
%         BW_Template = Canny_v0(Template_filt, filter, expe, 0, minIdx);
%         BW_Target = Canny_v0(Target_filt, filter, expe, 0, minIdx);
        BW_Template = edge(util_mask, filtro, [], sigmaValue);
        BW_Target = edge(geo_mask, filtro, [], sigmaValue);

%         if exist('filter','var')
%             if filter == 0
%                 BW_Template = edge(Template_filt, filtro, kValue*threshOut1, sigmaValue);
%                 BW_Target = edge(Target_filt, filtro, kValue*threshOut2, sigmaValue);
% 
%             elseif filter == 1
%                 BW_Template = Canny_v0(Template_filt, 1, kValue, expe, minIdx);
%                 BW_Target = Canny_v0(Target_filt, 1, kValue, expe, minIdx);
% 
%             elseif filter == 2
%                 BW_Template = Canny_v0(Template_filt, 2, kValue, expe, minIdx);
%                 BW_Target = Canny_v0(Target_filt, 2, kValue, expe, minIdx);
% 
%             elseif filter == 3
%                 BW_Template = Canny_v0(Template_filt, 3, kValue, expe, minIdx);
%                 BW_Target = Canny_v0(Target_filt, 3, kValue, expe, minIdx);
%                 
%             elseif filter == 4
%                 BW_Template = Canny_v0(Template_filt, 4, kValue, expe, minIdx);
%                 BW_Target = Canny_v0(Target_filt, 4, kValue, expe, minIdx);
% 
%             end
%         else
%             BW_Template = edge(Template_filt, filtro, kValue*threshOut1, sigmaValue,expe);
%             BW_Target = edge(Target_filt, filtro, kValue*threshOut2, sigmaValue,expe);
%         end
        
%         figure; montage({BW_Template, BW_Template1, BW_Template2, BW_Template3});
%         figure; subplot(221); imshow(BW_Target); subplot(222); imshow(BW_Target1);
%         subplot(223); imshow(BW_Target2); subplot(224); imshow(BW_Target3);
        
        dgb = 1;
        
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
%         k = 0.70;

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
    %         k = 0.70;

            BW_Template = edge(Template_filt, filtro, kValue(j)*threshOut1, sigmaValue);
            BW_Target = edge(Target_filt, filtro, kValue(j)*threshOut2, sigmaValue);
        end
        
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

% hFig = figure(55);
% hAx  = axes;
% imshow(Target,'Parent', hAx);
% imrect(hAx, [xoffSet, yoffSet, size(Template,2), size(Template,1)]);

centro = [ptCenterX ptCenterY];
end
