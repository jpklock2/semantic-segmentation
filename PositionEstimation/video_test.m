clear all;
close all;
clc;

%format long
arquivo = fopen('D:\payo\Documents\video4_SP_elcio\log_2020.txt');
temp = textscan(arquivo, '%s %s %s %s %s %s %s', [Inf Inf]);
fclose(arquivo); 
name=temp{1}; % nombre a
lat=str2double(temp{2}); %b
lon=str2double(temp{3}); %c
alt= str2double(temp{4}); %g
yaw=str2double(temp{5}); %d
pitch=str2double(temp{6});%e
roll=str2double(temp{7}); %f

tamx_calc = 400;
tamy_calc = 400;

log_coords = [];

[geo_img, cmap, R, bbox] = geotiffread('D:\payo\Documents\video4_SP_elcio\Mosaico.tif');

R_copy = R;

for i=1:length(name)
    
    filename1=('D:\payo\Documents\video4_SP_elcio\100CANON\');
    filename=sprintf('IMG_%d.JPG',(i+4562));  
    fprintf(filename);
    outf=fullfile(filename1,filename);
    uav_img = imread(outf);
    
    scale_uav_img = imresize(uav_img, 0.13, 'Bilinear'); %0.21
    figure(1);
    imshow(scale_uav_img);
    
    uav_size=size(scale_uav_img);
    uav_row_y=floor(uav_size(1)/2);
    uav_col_x=floor(uav_size(2)/2);
                                 
    hab=0;
    [crop_geo_img,cmap,R,bbox]=get_geo_subimg(lat(i,1),lon(i,1), ...
                                               tamx_calc,tamy_calc, ...
                                               geo_img, cmap, R_copy,... 
                                               bbox,hab);
    
    figure (2);
    imshow(crop_geo_img);

    [pre_geo_img, pre_uav_img]=preprocessing(crop_geo_img,scale_uav_img);
    
    [distx, disty, rotate_img]=correcion_perspectiva(...
       pre_uav_img, yaw(i), pitch(i), roll(i), alt(i));
    
    figure (3);
    imshowpair(crop_geo_img,rotate_img,'montage');
    figure (33);
    imshow(rotate_img);
    
    only_yaw_rotate_img = imrotate(pre_uav_img, -yaw(i));
    figure (4);
    imshow(only_yaw_rotate_img);
    
    
    [yoffSet,xoffSet,Mcorr,centro] = edges_and_correlation(...
                                      rotate_img, crop_geo_img);
    
    centro
    retorno = Mcorr;
    
         

    
    % Calculo da latitude e longitude com base nos pixeis da img.
    [lat_srp, lon_srp] = pix2latlon(R, centro(1), centro(2))
   
   
%     lat_p=sprintf('%0.10f', lat_srp)
%     lon_p=sprintf('%0.10f', lon_srp)
    
    dist = m_idist(lat_srp, ...
        lon_srp, ...
        lat(i,1), ...
        lon(i,1))
   % pause
    log_coords = [log_coords; lat(i,1),lon(i,1), lat_srp, lon_srp, dist];
    pause
end

FP = find(log_coords(:,5) > 100)

media=sum(log_coords(:,5))/length(log_coords(:,5))
