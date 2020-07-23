% IMAGENS DIFERENTE SENSOR
function [pre_geo_img, pre_uav_img]=preprocessing(geo_img,uav_img)

geo_img_copy = geo_img;
vant_img_copy = uav_img;

% figure(10);  imshowpair(geo_img_copy,vant_img_copy,'montage','Scaling',...
% 'joint');

img_vant_match = imhistmatch(vant_img_copy,geo_img_copy,256); 

if size(img_vant_match,3)==3
    img_vant_match=rgb2gray(img_vant_match);
else
    img_vant_match=vant_img_copy;
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

pre_geo_img = medfilt2(img_geo_eq_bl ,[3 3]);
pre_uav_img = medfilt2(img_vant_match_eq_bl_gray,[3 3]);
% figure(20);  imshowpair(pre_geo_img,pre_uav_img,'montage','Scaling','joint');
end