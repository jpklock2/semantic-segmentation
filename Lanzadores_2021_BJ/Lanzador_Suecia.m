%% Starting Code
initCode;
dataset = 'SC';
loadImageNames_BJ;
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
load('Images/Data_Suecia/frame_data.mat');
frameLat=FRAMEDATA.e01; % a lat
frameLong=FRAMEDATA.e1; %b lon
frameYaw=FRAMEDATA.e2; % c yaw
framePitch=FRAMEDATA.e00; %d pitch
frameRoll=FRAMEDATA.e3; % e roll
framePress=FRAMEDATA.e4; % f presion
%cargo los valores de posicion en pixels
load('Images/Data_Suecia/posicao_correta.mat');

a1=POSICAOCORRETA.e03; % lat
b1=POSICAOCORRETA.e1; %lon

load('Images/Data_Suecia/posicao_piobj.mat');
a3=POSICAOPIOBJ.e01;
b3=POSICAOPIOBJ.e00;


%% Reading Georeferenced Image
fprintf('\nProcessing Georeferenced Image...\n');
%Define el tamanho del recorte de la imagen georeferenciada 
tamx_calc = 600;
tamy_calc = 600;
[geo_img, cmap] = imread('Images/Data_Suecia/ortho_19comp.jpg');

%% Calculating Sensor Informations
%ref: https://www.imaging-resource.com/PRODS/T3I/T3IDAT.HTM
conv= pi/180;
fx=472.8380;
fy=518.0960;
ox=360/2;
oy=288/2;
camUpsideDown = 0;
imageRes=0.5;
px=(fx*imageRes)./framePress;
py=(fy*imageRes)./framePress;
imageRefLat=55.7238*conv;
imageRefLong=13.4712*conv;
imageRefX=1130.5;
imageRefY=936.5;
REARTH = 6400000;
incerteza_pixel = 80/imageRes;

log_coords1 = [];
templates = [{}];
targets = [{}];
allTemplateResults = [{}];

for i = 1:length(imageNamesTemp)-1
    %% Main Loop
    fprintf('\nRunning Image %d\n', i);
    
    %% Reading and Resizing Image
    fprintf('Reading and Resizing Image');
    filename1 = imagesPath;
    filename=sprintf('img%d.pgm',i); 
    outf = fullfile(filename1, filename);
    uav_img = imread(outf);
    
    [linhas,colunas,~]=size(uav_img); 
    dimX = linhas/px(i);
    dimY = colunas/py(i);
    scale_uav_img = imresize(uav_img, [dimX dimY], 'Bilinear');
    [linhas1,colunas1,~] = size(scale_uav_img);
    if printFigs
        figure;
        imshow(scale_uav_img);
    end
    
    %% Geo Sub Image
    fprintf('\nGetting Geo SubImage...\n');
    tic;
    dimWinEast = linhas1 + incerteza_pixel;
    dimWinNorth = colunas1 + incerteza_pixel;
    x_corte_ortho = (a1(i) - (dimWinEast/2));
    y_corte_ortho = (b1(i) - (dimWinNorth/2));
     
    crop_geo_img=imcrop(geo_img,[x_corte_ortho y_corte_ortho dimWinEast-1 dimWinNorth-1]);
    fprintf('Execution time for Geo SubImage: %f s\n', toc);
    
    %para efeitos de visualizaçao
    [col_geo_crop, lin_geo_crop,~] = size(crop_geo_img);
    pos= [col_geo_crop/2, lin_geo_crop/2, 3];
    color={'red'};
    crop_geo_img_withPoint= insertShape(crop_geo_img,'circle', pos, 'Color', color, 'LineWidth',3, 'Opacity',0.7);
    
    if printFigs
        figure;
        imshow(crop_geo_img);
        figure;
        imshow(crop_geo_img_withPoint);
    end
    
    %% Preprocessing
    fprintf('\nPreprocessing Image...\n');
    tic; 
    [pre_geo_img, pre_uav_img]=preprocessing_BJ(crop_geo_img,scale_uav_img,...
        filter_type_preprocessing, printFigs);   
    fprintf('Execution time for Preprocessing: %f s\n', toc);
    
    if printFigsCorr
        figure;
        imshowpair(pre_geo_img,pre_uav_img,'montage');
    end
    

    %% Perspective correction 
    fprintf('\nPerspective correction Image...\n');
    tic;
    [rot_img, util_rot_img, utilCropSize]=correcion_perspectiva_BJ(...
       scale_uav_img, (frameYaw(i)), framePitch(i), frameRoll(i), framePress(i),...
       pre_geo_img, pre_uav_img, perspective_correction_type, printFigs, preprocessing_img);
    fprintf('Execution time for Perspective Correction: %f s\n', toc);
    
    if printFigsCorr
        figure;
        imshowpair(crop_geo_img,util_rot_img,'montage');
    end
    
     %% Aplicando filtro
    if filter_img
        tic;
        rgbImage = util_rot_img;
        applyFilter_BJ;
        util_rot_img = rgbImage;
        clear rgbImage
        rgbImage = crop_geo_img;
        applyFilter_BJ;
        crop_geo_img = rgbImage;
        clear rgbImage
        rgbImage = pre_geo_img;
        applyFilter_BJ;
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
            
            [yoffSet,xoffSet,Mcorr,centro] = edges_and_correlation_BJ(...
                                      util_rot_img, pre_geo_img, n_case,printFigsCorr);
        else
            [yoffSet,xoffSet,Mcorr,centro] = edges_and_correlation_BJ(...
                                      util_rot_img, crop_geo_img, n_case,printFigsCorr);
        end
    
        fprintf('Execution time for Edges and Correlation: %f s\n', toc);

        if printFigsCorr
%             figure;
%             mesh(Mcorr);
        
            hFig = figure;
            hAx  = axes;
            crop_geo_img_withPoint1= insertShape(crop_geo_img_withPoint,'circle', [centro(1) centro(2) 3], 'Color', 'blue', 'LineWidth',3, 'Opacity',0.7);
            imshow(crop_geo_img_withPoint1,'Parent', hAx);
            imrect(hAx, [xoffSet, yoffSet, size(util_rot_img,2), size(util_rot_img,1)]);
        end
        
        match_x= a1(i)-dimWinEast/2 + centro(1);
        match_y= b1(i)-dimWinNorth/2 + centro(2);
        dX_match = (match_x-imageRefX)*imageRes;
        dY_match = -(match_y-imageRefY)*imageRes;
        Lat_match = (dY_match-a3(i))/REARTH + imageRefLat;
        Long_match = (dX_match-b3(i))/(REARTH*cos(imageRefLat))+ imageRefLong;
        erro_lat = (Lat_match-(frameLat(i)*conv))*REARTH;
        erro_long = (Long_match-(frameLong(i)*conv))*REARTH*cos(imageRefLat);
        ERR_TOTAL = sqrt(power(erro_lat,2) + power(erro_long,2))
       

        log_coords1 = [log_coords1; ERR_TOTAL,n_case,i];       
%         save log_coords1.mat log_coords1

    end    
end

% FP = find(log_coords1(:,5) > 100)
% 
% 
% media=sum(log_coords1(:,5))/length(log_coords1(:,5));

case1 = log_coords1(find(log_coords1(:,2) == 1),1);
mean_case1 = mean(case1)

case2 = log_coords1(find(log_coords1(:,2) == 2),1);
mean_case2 = mean(case2)

case3 = log_coords1(find(log_coords1(:,2) == 3),1);
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