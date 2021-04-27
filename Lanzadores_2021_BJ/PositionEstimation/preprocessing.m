% IMAGENS DIFERENTE SENSOR
function [pre_geo_img, pre_uav_img]=preprocessing(geo_img,uav_img,...
    filter_type, printFigs)

geo_img_copy = geo_img;
uav_img_copy = uav_img;

if size(geo_img_copy,3)~=3
    geo_img_copy = cat(3, geo_img_copy, geo_img_copy, geo_img_copy);
end
if size(uav_img_copy,3)~=3
    uav_img_copy = cat(3, uav_img_copy, uav_img_copy, uav_img_copy);
end

img_vant_match = imhistmatch(uav_img_copy,geo_img_copy,256);   
if size(img_vant_match,3)==3
    img_vant_match=rgb2gray(img_vant_match);
% else
%     img_vant_match=img_vant_match;
end

img_vant_match_eq = histeq(img_vant_match);
if size(geo_img_copy,3)==3
    img_geo_gray=rgb2gray(geo_img_copy);
else
    img_geo_gray=geo_img_copy;
end

img_geo_eq = histeq(img_geo_gray);
img_vant_match_eq_gray = imhistmatch(img_vant_match_eq,img_geo_eq,256);

w=fspecial('laplacian',0);
f=im2double(img_geo_eq);
gw=imfilter(f,w,'replicate');
img_2=f-gw;
img_geo_eq_bl = medfilt2(img_2,[2 2]);

w=fspecial('laplacian',0);
f_1=im2double(img_vant_match_eq_gray);
gw_1=imfilter(f_1,w,'replicate');
img_1=f_1-gw_1;
img_vant_match_eq_bl_gray = medfilt2(img_1,[2 2]);

img_vant_match_eq_bl_gray = imhistmatch(img_vant_match_eq_bl_gray,...
    img_geo_eq_bl,256);

if filter_type == 1 % filter_type = 'media'
    pre_geo_img = medfilt2(img_geo_eq_bl ,[3 3]); % Iin
    pre_uav_img = medfilt2(img_vant_match_eq_bl_gray,[3 3]); %Iout
    
    if printFigs
        figure;
        imshowpair(pre_geo_img,pre_uav_img,'montage','Scaling','joint');
    end
    
elseif filter_type == 2 % filter_type = 'ILS'
    lambda = 0.6;
    iter = 3;
    p = 0.3;
    eps = 0.0001;  

    Img_geo_ILS = im2double(img_geo_eq_bl);
    Smoothed_geo = ILS_LNorm(Img_geo_ILS, lambda, p, eps, iter);
    Diff_geo = Img_geo_ILS - Smoothed_geo;
    pre_geo_img = Img_geo_ILS + 3 * Diff_geo;

    Img_uav_ILS = im2double(img_vant_match_eq_bl_gray);
    Smoothed_uav = ILS_LNorm(Img_uav_ILS, lambda, p, eps, iter);
    Diff_uav = Img_uav_ILS - Smoothed_uav;
    pre_uav_img = Img_uav_ILS + 3 * Diff_uav;
    
    if printFigs
        figure;
        imshowpair(pre_geo_img,pre_uav_img,'montage','Scaling','joint');
    end
end
pre_geo_img = im2uint8(pre_geo_img);
pre_uav_img = im2uint8(pre_uav_img);

end