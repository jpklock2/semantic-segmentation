%% Starting Code
initCode;
imageNamesTemp = imageNames(2:end);

%% Reading Logs
arquivo = fopen('Images\logs.txt');
temp = textscan(arquivo, '%s %s %s %s %s %s %s', [Inf Inf]);
fclose(arquivo); 
name=temp{1}; % nombre a
idx = [];
for i = 1:length(imageNamesTemp)
    idx = [idx; find(contains(name,imageNamesTemp{i}))];
end
lat=str2double(temp{2}(idx)); %b
lon=str2double(temp{3}(idx)); %c
alt= str2double(temp{4}(idx)); %g
yaw=str2double(temp{5}(idx)); %d
pitch=str2double(temp{6}(idx));%e
roll=str2double(temp{7}(idx)); %f

%% Reading Georeferenced Image
fprintf('\nProcessing Georeferenced Image...\n');
tamx_calc = 400;
tamy_calc = 400;
log_coords = [];
[geo_img, cmap, R, bbox] = geotiffread('Images\Train\Mosaicoo.tif');
R_copy = R;

%% Training Georeferenced Image
fprintf('\nRunning Semantic Segmentation Pipeline...\n');
fprintf('\n%s\n', repmat('º', [1 75]));
[maskGeo, maskIdxGeo, centroidsGeo, classesGeo, parameters] = mainFunc('Images\Train\Mosaicoo.tif', 1);
fprintf('%s\n', repmat('º', [1 75]));
fprintf('\n%s\n', repmat('-', [1 100]));

printFigs = 0;
for i=1:length(imageNamesTemp)
    %% Main Loop
    fprintf('\nRunning Image %d\n', i);
    
    %% Reading and Resizing Image
    filename1=('Images\Test\Original');
    outf=fullfile(filename1,strtrim(imageNamesTemp{i}));
    uav_img = imread(outf);
    scale_uav_img = imresize(uav_img, 0.13, 'Bilinear'); %0.21
    if printFigs
        figure(1);
        imshow(scale_uav_img);
    end
    uav_size=size(scale_uav_img);
    uav_row_y=floor(uav_size(1)/2);
    uav_col_x=floor(uav_size(2)/2);
    
    %% Geo Sub Image
    fprintf('\nGetting Geo SubImage...\n');
    tic;
    hab=0;
    [crop_geo_img,cmap,R,bbox]=get_geo_subimg(lat(i,1),lon(i,1), ...
                                               tamx_calc,tamy_calc, ...
                                               geo_img, cmap, R_copy,... 
                                               bbox,hab);
    if printFigs                                       
        figure (2);
        imshow(crop_geo_img);
    end
    fprintf('Execution time for Geo SubImage: %f s\n', toc);
    
    %% Semantic Mask
    fprintf('\nRunning Semantic Segmentation Pipeline...\n');
    fprintf('\n%s\n', repmat('º', [1 75]));
    [mask, maskIdx, centroids, classes] = mainFunc(outf, i+1, classesGeo, centroidsGeo, parameters);
    fprintf('%s\n', repmat('º', [1 75]));

    %% Preprocessing
    fprintf('\nPreprocessing Image...\n');
    tic;
    [pre_geo_img, pre_uav_img]=preprocessing(crop_geo_img,scale_uav_img);
    fprintf('Execution time for Preprocessing: %f s\n', toc);
    
    %% Correcting Perspective
    fprintf('\nCorrecting Image Perspective...\n');
    tic;
    [distx, disty, rotate_img]=correcion_perspectiva(...
       pre_uav_img, yaw(i), pitch(i), roll(i), alt(i));
    
    if printFigs
        figure (3);
        imshowpair(crop_geo_img,rotate_img,'montage');
        figure (33);
        imshow(rotate_img);
    end
    only_yaw_rotate_img = imrotate(pre_uav_img, -yaw(i));
    if printFigs
        figure (4);
        imshow(only_yaw_rotate_img);
    end
    fprintf('Execution time for Correcting Perspective: %f s\n', toc);
    
    %% Edges and Correlation
    fprintf('\nGetting Edges and Correlation...\n');
    tic;
    [yoffSet,xoffSet,Mcorr,centro] = edges_and_correlation(...
                                      rotate_img, crop_geo_img);
    
    centro;
    retorno = Mcorr;
    fprintf('Execution time for Edges and Correlation: %f s\n', toc);
    
    %% Distance
    fprintf('\nCalculating Distance...\n');
    tic;
    % Calculo da latitude e longitude com base nos pixeis da img.
    [lat_srp, lon_srp] = pix2latlon(R, centro(1), centro(2));
    
%     lat_p=sprintf('%0.10f', lat_srp)
%     lon_p=sprintf('%0.10f', lon_srp)
    
    dist = m_idist(lat_srp, ...
        lon_srp, ...
        lat(i,1), ...
        lon(i,1));
    fprintf('Execution time for Distance: %f s\n', toc);
    
    %% Saving Results
    log_coords = [log_coords; lat(i,1),lon(i,1), lat_srp, lon_srp, dist];
%     pause
    fprintf('\n%s\n', repmat('-', [1 100]));
end

%% Final Results
FP = find(log_coords(:,5) > 100);
media=sum(log_coords(:,5))/length(log_coords(:,5));

fprintf('\n');
