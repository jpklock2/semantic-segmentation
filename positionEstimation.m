%% Starting Code
initCode;
useAll = 0;
baseCase = 0;
loadImageNames;
imageNamesTemp = imageNames(2:end);
startTime = datetime('now');
fprintf('\nExperiment start: %s\n', startTime);

%% Reading Logs
arquivo = fopen('Images\logs.txt');
temp = textscan(arquivo, '%s %s %s %s %s %s %s', [Inf Inf]);
fclose(arquivo);
name=temp{1}; % nombre a
idx = [];
for i = 1:length(imageNamesTemp)
    for j = 1:length(name)
        if contains(imageNamesTemp{i}, name{j})
            idx = [idx; j];
            break;
        end
    end
end
lat=str2double(temp{2}(idx)); %b
lon=str2double(temp{3}(idx)); %c
alt= str2double(temp{4}(idx)); %g
yaw=str2double(temp{5}(idx)); %d
pitch=str2double(temp{6}(idx));%e
roll=str2double(temp{7}(idx)); %f

%% Reading Georeferenced Image
fprintf('\nProcessing Georeferenced Image...\n');
tamx_calc = 800;
tamy_calc = 800;
log_coords1_all = [{}];
log_coords2_all = [{}];
[geo_img, cmap, R, bbox] = geotiffread('Images\Train\Mosaicoo.tif');
R_copy = R;

% labImage = rgb2lab(geo_img);
% L = labImage(:,:,1)/100;
% L = adapthisteq(L);
% labImage(:,:,1) = L*100;
% geo_img = lab2rgb(labImage);
% 
% geo_img_rot = imrotate(geo_img, 4.04);
% geo_img_crop = imcrop(geo_img_rot, [5913.5 2679.5 1083 722]);
% imwrite(geo_img_crop, ['E:/semantic-segmentation/Images/Test/Cropped/cropped_image_4.JPG']);
% figure; imshow(geo_img_crop);

%% Calculating Sensor Informations
%ref: https://www.imaging-resource.com/PRODS/T3I/T3IDAT.HTM
focal_length_mm = 35;
sensor_w_mm = 22.3;
sensor_h_mm = 14.9;
size_w_uav_img = 5184;
size_h_uav_img = 3456;

focal_length_px = (focal_length_mm*size_w_uav_img)/sensor_w_mm;
focal_length_py = (focal_length_mm*size_h_uav_img)/sensor_h_mm;
ox = size_w_uav_img/2;
oy = size_h_uav_img/2;
imgRes = 0.5;

log_coords1_all = [{}];
log_coords2_all = [{}];
templates = [{}];
targets = [{}];

allTemplateResults = [{}];
% load('templateResults.mat')

% expe = 0.2:0.2:2;

% load('experiment_temp.mat');
% size(log_coords1_all, 1)+
printFigs = 0;
for i = size(allTemplateResults, 1)+1:length(imageNamesTemp)
    
%     for expe = 0.2:0.2:2
    for expe = 1:1
        
    %% Main Loop
    fprintf('\nRunning Image %d\n', i);
    
    %% Reading and Resizing Image
    filename1 = imagesPath;
    outf = fullfile(filename1,strtrim(imageNamesTemp{i}));
    uav_img = imread(outf);
    
    getHomData = 2;
    if getHomData == 0
        [movingPoints,fixedPoints] = cpselect(uav_img(end:-1 : 1, end:-1 : 1, :),geo_img,'Wait',true);
        tform = fitgeotrans(movingPoints,fixedPoints,'projective');
        load('transformData.mat');
        dimX = size(uav_img, 1)*tform.T(1,1);
        dimY = size(uav_img, 2)*tform.T(2,2);
    elseif getHomData == 1
%         scale_uav_img = imresize(uav_img, 0.13, 'Bilinear'); %0.21
        dimX = size(uav_img, 1)*0.206944444444444;
        dimY = size(uav_img, 2)*0.206944444444444;
    elseif getHomData == 2
        [linhas,colunas,~]=size(uav_img); 
        px=(focal_length_px*imgRes)./alt(i,1);
        py=(focal_length_py*imgRes)./alt(i,1);
        dimX = linhas/px;
        dimY = colunas/py;
    end
    
    if useCropped
        outf2 = fullfile('Images/Test/Original_Dev',strtrim(imageNamesTemp{i}(end-11:end)));
        uav_img_temp = imread(outf2);
        [linhas,colunas,~]=size(uav_img_temp); 
        dimX = linhas/px;
        dimY = colunas/py;
    end
    scale_uav_img = imresize(uav_img, [dimX dimY], 'Bilinear');
    if printFigs
        figure(1);
        imshow(scale_uav_img);
    end
    
    %% Semantic Mask
    
    %     cases = [1 2 5];
    cases = 2;
    n_case = 1;
%     filters = 3:5;
    filters = 3;
    
    log_coords1 = [];
    log_coords2 = [];
    
%     for n_case = 1:length(cases)
    for n_filt = 1:length(filters)
    
    if ~baseCase
        
    %% Training Georeferenced Image
    fprintf('\nRunning Semantic Segmentation Pipeline...\n');
    fprintf('\n%s\n', repmat('º', [1 75]));
    processGeo = 0;
    if exist(['geoData_filter_' num2str(filters(n_filt)) '.mat'], 'file') == 2 && ~processGeo
        fprintf('\nGeoreferenced data found, loading data...\n\n');
        load(['geoData_filter_' num2str(filters(n_filt)) '_25.mat']);
    else
        [maskGeo, maskIdxGeo, centroidsGeo, classesGeo, parameters, geoAdjacencies, ftGeoOwn, ftGeoAdj] = semanticSegmentation('Images\Train\Mosaicoo.tif', 1, filters(n_filt));
        LmaskGeo = zeros(size(maskGeo));
        for mx = 1:length(classesGeo)
            LmaskGeo(maskGeo == mx) = classesGeo(mx);
        end
        save(['geoData_filter_' num2str(filters(n_filt)) '_25.mat'],'maskGeo','maskIdxGeo','centroidsGeo','classesGeo','parameters','LmaskGeo','geoAdjacencies','ftGeoOwn','ftGeoAdj');
    end

%     end
%     end
%     for n_filt = 1:length(filters)

%     colorLmask = zeros(size(maskGeo));
%     for cl = 1:length(maskIdxGeo)
%         colorLmask(maskIdxGeo{cl}) = classesGeo(cl);
%     end
%     figure; imagesc(colorLmask);
% 
%     colorSuperpixelFunc(geo_img, classesGeo, maskIdxGeo, parameters);
%     textureAnalysis(parameters, centroidsGeo, maskGeo, maskIdxGeo, classesGeo, geo_img);

    fprintf('%s\n', repmat('º', [1 75]));

    fprintf('\n%s\n', repmat('-', [1 100]));
    
%         figure; imagesc(Lmask);
        
        fprintf('\nRunning Semantic Segmentation Pipeline...\n');
        fprintf('\n%s\n', repmat('º', [1 75]));
        [mask, maskIdx, centroids, classes, ~, adjacencies, ftOwn, ftAdj, currImage, currImagePlot] = semanticSegmentation(outf, i+1, filters(n_filt), classesGeo, centroidsGeo, parameters, [dimX dimY]);
        currImage = imresize(currImage, [dimX dimY], 'Bilinear');
        outputSegmentation = evalFunction(classes, length(unique(classes)), maskIdx, currImagePlot, length(classes));
        %         oldMask = mask;
        Lmask = zeros(size(mask));
        for mx = 1:length(classes)
            Lmask(mask == mx) = classes(mx);
        end
%         load(['geoData_filter_' num2str(filters(n_filt)) '.mat']);
%         figure; imagesc(Lmask);
        fprintf('%s\n', repmat('º', [1 75]));
    end
    
    %% Geo Sub Image
    fprintf('\nGetting Geo SubImage...\n');
    tic;
    hab=1;
    [crop_geo_img,cmap,R,bbox,cropSize]=get_geo_subimg(lat(i,1),lon(i,1),...
        tamx_calc, tamy_calc, geo_img, cmap, R_copy, bbox, hab);
    
    if printFigs                                       
        figure (2);
        imshow(crop_geo_img);
    end
    fprintf('Execution time for Geo SubImage: %f s\n', toc);

    %% Preprocessing
    fprintf('\nPreprocessing Image...\n');
    tic;
    [pre_geo_img, pre_uav_img]=preprocessing(crop_geo_img,scale_uav_img);
    fprintf('Execution time for Preprocessing: %f s\n', toc);
    
    %% Correcting Perspective
    fprintf('\nCorrecting Image Perspective...\n');
    tic;
    if ~baseCase
        mask = imresize(mask, [size(pre_uav_img, 1) size(pre_uav_img, 2)]);
        Lmask = imresize(Lmask, [size(pre_uav_img, 1) size(pre_uav_img, 2)]);
        [rot_img, util_rot_img, rot_mask, util_rot_mask, util_mask, utilCropSize, currImage, rotPlotImage, rotOutputSegmentation]=correcion_perspectiva(...
            scale_uav_img, (yaw(i)-20), pitch(i), roll(i), alt(i), Lmask, mask, currImage, currImagePlot, outputSegmentation);
    else
        [rot_img, util_rot_img]=correcion_perspectiva(...
            scale_uav_img, (yaw(i)-20), pitch(i), roll(i), alt(i));
    end
   
    if ~baseCase
        % removing residuals from resize
        rot_mask = removeResiduals(rot_mask, classes);
        util_rot_mask = removeResiduals(util_rot_mask, classes);
        util_mask = removeResidualsMask(util_mask);

%         figure; montage(rot_img, imagesc(rot_mask));
    end
    
    if printFigs
        figure (3);
        imshowpair(crop_geo_img,rot_img,'montage');
    end

    if printFigs
        figure (4);
        imshowpair(crop_geo_img,util_rot_img,'montage');
    end
    fprintf('Execution time for Correcting Perspective: %f s\n', toc);
    
    %% Get Useful Area
    
%     only_yaw_rotate_img = imrotate(scale_uav_img, -yaw(i));
%     uav_size = size(scale_uav_img);
%     util_only_yaw_rotate_img = useful_area(...
%         only_yaw_rotate_img, uav_size(1), uav_size(2),-yaw(i));
    
    %% Edges and Correlation
    
    run_edges_and_correlation = 0;
    if run_edges_and_correlation
    fprintf('\nGetting Edges and Correlation...\n');
    fprintf('Case %d:\n', cases(n_case));
    
    tic;
    if ~baseCase
        [yoffSet,xoffSet,Mcorr,centro,kVal,kVal2,kVal3,kVal4,bwTemplate,bwTarget] = edges_and_correlation(...
            util_rot_img, crop_geo_img, cases(n_case), rot_mask, util_rot_mask, util_mask, Lmask, LmaskGeo, centroids, classes, classesGeo,...
            parameters, filters(n_filt), expe, crop_geo_mask, adjacencies, ftOwn, ftAdj, geoAdjacencies, ftGeoOwn,...
            ftGeoAdj, maskGeo, cropSize, utilCropSize, currImage, geo_img, centroidsGeo);
    else
        [yoffSet,xoffSet,Mcorr,centro] = edges_and_correlation(...
            util_rot_img, crop_geo_img, cases(n_case));
    end
    centro;
    retorno = Mcorr;
    templates = [templates; bwTemplate];
    targets = [targets; bwTarget];
    fprintf('Execution time for Edges and Correlation: %f s\n', toc);
    
    %% Distance
    fprintf('\nCalculating Distance...\n');
    tic;
    % Calculo da latitude e longitude com base nos pixeis da img.
    [lat_srp, lon_srp] = pix2latlon(R, centro(1), centro(2));
    [lat_srp2, lon_srp2] = pix2latlon(R, centro(2), centro(1));
    
%     lat_p=sprintf('%0.10f', lat_srp)
%     lon_p=sprintf('%0.10f', lon_srp)
    
    dist1 = m_idist(lon_srp, lat_srp, lon(i,1), lat(i,1));
    dist2 = m_idist(lon_srp2, lat_srp2, lon(i,1), lat(i,1));
    
    fprintf('Execution time for Distance: %f s\n', toc);
    
    %% Saving Results
%     log_coords1 = [log_coords1; lat(i,1),lon(i,1), lat_srp, lon_srp, dist1, cases(n_case), i];
%     log_coords2 = [log_coords2; lat(i,1),lon(i,1), lat_srp2, lon_srp2, dist2, cases(n_case), i];
    log_coords1 = [log_coords1; lat(i,1),lon(i,1), lat_srp, lon_srp, dist1, filters(n_filt), kVal, kVal2, kVal3, kVal4, i];
    log_coords2 = [log_coords2; lat(i,1),lon(i,1), lat_srp2, lon_srp2, dist2, filters(n_filt), kVal, kVal2, kVal3, kVal4, i];
    
    else
        [resCentro, resProb, resIn, resArea, resDist, visDist, visProb, visIn, redArea] = geoAdjacencyTM(...
            geoAdjacencies, adjacencies, classes, Lmask, util_mask, LmaskGeo,...
            cropSize, utilCropSize, parameters, filters(n_filt), geo_img, util_rot_mask, R_copy, lon(i,1), lat(i,1), i, rotPlotImage, rotOutputSegmentation);
        
        [yoffSet, xoffSet, Mcorr, centro] = edges_and_correlation_v2(util_rot_img, crop_geo_img, cases(n_case), resArea, cropSize, redArea);
        
        centro;
        retorno = Mcorr;
%         templates = [templates; bwTemplate];
%         targets = [targets; bwTarget];
        
        %% Distance
        fprintf('\nCalculating Distance...\n');
        tic;
        % Calculo da latitude e longitude com base nos pixeis da img.
        [lat_srp, lon_srp] = pix2latlon(R, centro(1), centro(2));
        [lat_srp2, lon_srp2] = pix2latlon(R, centro(2), centro(1));

    %     lat_p=sprintf('%0.10f', lat_srp)
    %     lon_p=sprintf('%0.10f', lon_srp)

        dist1 = m_idist(lon_srp, lat_srp, lon(i,1), lat(i,1));
        dist2 = m_idist(lon_srp2, lat_srp2, lon(i,1), lat(i,1));

        fprintf('Execution time for Distance: %f s\n', toc);
        
    end
    
    end
    
%     log_coords1_all = [log_coords1_all; log_coords1];
%     log_coords2_all = [log_coords2_all; log_coords2];
%     save('experiment_temp.mat', 'log_coords1_all', 'log_coords2_all');
%     pause

    fprintf('\n%s\n', repmat('-', [1 100]));
    
    end
    
    allTemplateResults = [allTemplateResults; {resCentro, resProb, resIn, resArea, resDist, visDist, visProb, visIn, centro, dist1}];
    save templateResults.mat allTemplateResults

% end

endTime = datetime('now');
fprintf('\nExperiment end: %s\n', endTime);
chartStart = replace(char(startTime), ':', '-');
charEnd = replace(char(endTime), ':', '-');
% save(['experiment_' chartStart(end-7:end) '_' charEnd(end-7:end) '.mat'], 'log_coords1_all', 'log_coords2_all');

end

%% Final Results
FP = find(log_coords1(:,5) > 100);
media=sum(log_coords1(:,5))/length(log_coords1(:,5));

% Métricas
log_coords = log_coords1_all;
allDists = [];
error_lat = [];
error_lon = [];
EPE = [];
EGPE = [];
Rearth = 6.3781*1e6;
for i = 1:size(log_coords, 1)
    allDists = [allDists log_coords{i}(:,5)];
    error_lat = [error_lat abs((log_coords{i}(:,3)-log_coords{i}(:,1)).*Rearth)];
    error_lon = [error_lon abs((log_coords{i}(:,4)-log_coords{i}(:,2)).*Rearth.*cos(log_coords{i}(:,1)))];
    EPE = [EPE sqrt(error_lat(:, end).^2 + error_lon(:, end).^2)];
    EGPE = [EGPE (log_coords{i}(:,5) <= 10)];
end
SEPE = [];
for i = 1:length(allDists)-29
    SEPE = [SEPE std(allDists(:, i:i+29), [], 2)];
end

% load('Temp\datasetResults.mat');

% cases = [1 2 7];
for i = 1:4
% testData = 100*rand([1 20]);
% testData = [allDists(:, 1); allDists(:, 2); allDists(:, 3)];
testData = allDists(i, :);
% testData = 100*rand([1 2000]);
% testData = [linspace(1, 10) linspace(1, 10) linspace(1, 10)];
x = 1:length(testData);

figure;
clear g
g = gramm('x',x,'y',testData,'color',testData); % cria estrutura
g.set_continuous_color('LCH_colormap',[20 80 ; 40 30 ; 260 260 ]); % barra de cores
g.geom_point('dodge', 0.8, 'alpha', 0.05); % cria pontos
g.stat_smooth(); % cria intervalo de confiança
% g.stat_smooth('method','loess'); % cria intervalo de confiança
g.set_color_options('map', 'matlab'); % muda espaço de cores
g.set_stat_options('alpha', 0.05); % tamanho do intervalo de confiança
% g.stat_summary('geom',{'bar'},'dodge',0,'width',0.4);
g.geom_bar('stacked', true, 'dodge', 0.8, 'FaceAlpha', 0.8, 'EdgeAlpha', 0); % cria barras
g.set_title(['Case ' num2str(cases(i)) ' results']); % titulo
g.draw(); % desenha
set([g.results.stat_smooth.area_handle],'FaceColor',[0.4 0.4 0.4]); % muda área
set([g.results.geom_bar_handle],'FaceColor','flat'); % muda cor das barras
set([g.results.geom_bar_handle],'FaceVertexCData',[0; 1]); % muda limite de cor das barras
set([g.results.geom_bar_handle],'CData',g.results.geom_point_handle.CData); % cor dos pontos pras barras
set([g.results.stat_smooth.line_handle],'Color','blue'); % cor da linha

figure;
[f,xi] = ksdensity(testData);
plot(xi,f, 'Linewidth', 2); hold on;
xiM = mean(xi);
fM = interp1(xi,f,xiM);
plot([xiM xiM], [0 fM], 'Linewidth', 2, 'Color', 'k', 'LineStyle', '--');
title(['Case ' num2str(cases(i)) ' density']);
legend('Density Curve', ['Mean = ' num2str(round(xiM, 2))]);
xlabel('EPE'); ylabel('Density');
% clear g
% g = gramm('x',x,'y',testData); % cria estrutura
% g.stat_density('function','pdf');
% g.set_color_options('map', 'matlab'); % muda espaço de cores
% g.set_title(['Case ' num2str(cases(i)) ' results']); % titulo
% g.draw(); % desenha

end

fprintf('\n');
