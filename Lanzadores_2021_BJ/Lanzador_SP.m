%% Starting Code
initCode;
dataset = 'SP';
loadImageNames;
imageNamesTemp = imageNames(2:end);
startTime = datetime('now');
fprintf('\nExperiment start: %s\n', startTime);

%% Configuration Code
printFigs = 0; %mostrar figuras
printFigsCorr = 0; %mostrar figuras de correlaçao
filter_img = 1; %aplica filtro sobre as imagens UAV e Geo
myFilter = 4; % 1 NO FILTER | 3 BILATERAL | 2 KUWAHARA | 5 ANISIOTROPIC | 4 ILS - Tipos de filtros para reduzir textura
preprocessing_img = 1; %se o template matching vai com imagem preprocessadas ou nao
filter_type_preprocessing = 2; % 1 MEDIA FILTER | 2 ILS - tipo de filtro
perspective_correction_type = 1; % 1 PARAMETRIC_LOG | 2 KEY_POINTS - tipo de corecçao de perspertiva
edge_type = 1; %futuramente para escolher o melhor extrator de borda

%% Reading Logs
fprintf('\nReading Logs... \n');
% arquivo = fopen('Images\Data_BR_SP\log_BR_SP.txt');
arquivo = fopen('Images/Data_BR_SP/log_BR_SP.txt'); %linux
temp = textscan(arquivo, '%s %s %s %s %s %s %s', [Inf Inf]);
fclose(arquivo); 
name=temp{1};
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
tamx_calc = 300; 
tamy_calc = 300;
log_coords1 = [];
log_coords2 = [];
% [geo_img, cmap, R, bbox] = geotiffread('Images\Data_BR_SP\Mosaico.tif');
[geo_img, cmap, R, bbox] = geotiffread('Images/Data_BR_SP/Mosaico.tif');
R_copy = R;

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
imgRes= 0.5;

log_coords1_all = [{}];
log_coords2_all = [{}];
templates = [{}];
targets = [{}];
allTemplateResults = [{}];

for i = 1:length(imageNamesTemp)
    
    %% Main Loop
    fprintf('\nRunning Image %d\n', i);
    
    %% Reading and Resizing Image
    filename1 = imagesPath;
    outf = fullfile(filename1, strtrim(imageNamesTemp{i}));
    uav_img = imread(outf);
    
    [linhas,colunas,~]=size(uav_img); 
    px=(focal_length_px*imgRes)./alt(i,1);
    py=(focal_length_py*imgRes)./alt(i,1);
    dimX = linhas/px;
    dimY = colunas/py;
    scale_uav_img = imresize(uav_img, [dimX dimY], 'Bilinear');
    scale_uav_img = imresize(scale_uav_img, 0.87, 'Bilinear');
    if printFigs
        figure;
        imshow(scale_uav_img);
    end
    
    %% Geo Sub Image
    fprintf('\nGetting Geo SubImage...\n');
    tic;
    hab=1;
    [crop_geo_img,crop_geo_img_withPoint,cmap,R,bbox,cropSize]=get_geo_subimg(lat(i,1),lon(i,1),...
        tamx_calc, tamy_calc, geo_img, cmap, R_copy, bbox, hab);
    
    fprintf('Execution time for Geo SubImage: %f s\n', toc);
    
    if printFigs
        figure;
        imshow(crop_geo_img);
    end
    
    %% Preprocessing
    fprintf('\nPreprocessing Image...\n');
    tic; 
    [pre_geo_img, pre_uav_img]=preprocessing(crop_geo_img,scale_uav_img,...
        filter_type_preprocessing, printFigs);   
    fprintf('Execution time for Preprocessing: %f s\n', toc);
    
    if printFigsCorr
        figure;
        imshowpair(pre_geo_img,pre_uav_img,'montage');
    end
    

    %% Correcting Perspective
    fprintf('\nCorrecting Image Perspective...\n');
%     tic;
%     
%     [rot_img, util_rot_img, utilCropSize]=correcion_perspectiva(...
%        scale_uav_img, (yaw(i)-20), pitch(i), roll(i), alt(i),...
%        pre_geo_img, pre_uav_img, perspective_correction_type, printFigs, preprocessing_img);
%     fprintf('Execution time for Perspective Correction: %f s\n', toc);
    
%     
%     rot_img = imrotate(scale_uav_img, -(yaw(i)-20),'Bilinear','crop');
%     [x2, y2] = size(rot_img);
%     util_rot_img =rot_img(abs(x2/2-195):abs(x2/2+195), abs(y2/2-195):abs(y2/2+195));
    
    [rot_img, util_rot_img, utilCropSize]=correcion_perspectiva_sem(...
        scale_uav_img, (yaw(i)-20), pitch(i), roll(i), alt(i),...
        pre_geo_img, pre_uav_img, perspective_correction_type, printFigs,...
        preprocessing_img);
    [x2, y2] = size(rot_img);
    util_rot_img =rot_img(abs(x2/2-195):abs(x2/2+195), abs(y2/2-195):abs(y2/2+195));

    
     %% Aplicando filtro
    if filter_img
        tic;
        rgbImage = util_rot_img;
        applyFilter;
        util_rot_img = rgbImage;
        clear rgbImage
        rgbImage = crop_geo_img;
        applyFilter;
        crop_geo_img = rgbImage;
        clear rgbImage
        rgbImage = pre_geo_img;
        applyFilter;
        pre_geo_img = rgbImage;
        clear rgbImage
        fprintf('\nPerspective correction Image: %f s\n', toc);
    end
    
    %% Edges and Correlation
    
    for n_case=1:3
        fprintf('\nGetting Edges and Correlation Matrix...\n');
        fprintf('Case %d:\n', n_case);
        tic;
        if preprocessing_img 
            
            [yoffSet,xoffSet,Mcorr,centro] = edges_and_correlation(...
                                      util_rot_img, pre_geo_img, n_case,printFigsCorr);
        else
            [yoffSet,xoffSet,Mcorr,centro] = edges_and_correlation(...
                                      util_rot_img, crop_geo_img, n_case,printFigsCorr);
        end
    
        fprintf('Execution time for Edges and Correlation: %f s\n', toc);

        if printFigsCorr
%             figure;
%             mesh(Mcorr);
        
            hFig = figure;
            hAx  = axes;
            crop_geo_img_withPoint1= insertShape(crop_geo_img_withPoint,'circle', [centro(1) centro(2) 30], 'Color', 'blue', 'LineWidth',25, 'Opacity',0.7);
            imshow(crop_geo_img_withPoint1,'Parent', hAx);
            imrect(hAx, [xoffSet, yoffSet, size(util_rot_img,2), size(util_rot_img,1)]);
        
        end

        % Calculo da latitude e longitude com base nos pixeis da img.
        [lat_srp, lon_srp] = pix2latlon(R, centro(1), centro(2));
 
        dist1 = m_idist(lon_srp, lat_srp, lon(i,1), lat(i,1))
        

        log_coords1 = [log_coords1; lat(i,1),lon(i,1), lat_srp, lon_srp,...
                       dist1,n_case,i];
        
%         save log_coords1.mat log_coords1

    end    
end

FP = find(log_coords1(:,5) > 100)


media=sum(log_coords1(:,5))/length(log_coords1(:,5));

case1 = log_coords1(find(log_coords1(:,6) == 1),5);
mean_case1 = mean(case1)

case2 = log_coords1(find(log_coords1(:,6) == 2),5);
mean_case2 = mean(case2)

case3 = log_coords1(find(log_coords1(:,6) == 3),5);
mean_case3 = mean(case3)
% 
% case4 = log_coords1(find(log_coords1(:,6) == 4),5);
% mean_case4 = mean(case4)

% case5 = log_coords1(find(log_coords1(:,6) == 5),5);
% mean_case5 = mean(case5)
% 
% case6 = log_coords1(find(log_coords1(:,6) == 6),5);
% mean_case6 = mean(case6)
% 
% case7 = log_coords1(find(log_coords1(:,6) == 7),5);
% mean_case7 = mean(case7)
% 
% case8 = log_coords1(find(log_coords1(:,6) == 8),5);
% mean_case8 = mean(case8)
% 
% case9 = log_coords1(find(log_coords1(:,6) == 9),5);
% mean_case9 = mean(case9)

figure; boxplot([case1,case2,case3])
