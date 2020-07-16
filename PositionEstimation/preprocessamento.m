close all;
clear;
clc

% IMAGENS DIFERENTE SENSOR
for ii=1:5
    fprintf('Img_teste: %d \n',ii);
    % Leitura das imagens
    filename1=('C:\Users\payo\Dropbox\COPY\siftDemoV4\Teste_paper_16_02_2018\correcciones_IEEE_ACCESS_28_09_2019\img_dif_sensor\');
    filename=sprintf('IMG_v_%d.png',ii);  
    outf=fullfile(filename1,filename);
    % Imagem do VANT
    img_vant = imread(outf);
    figure(1); imshow(img_vant)
    filename11=sprintf('IMG_g_%d.png',ii);  
    outf1=fullfile(filename1,filename11);
    % Imagem Georeferenciada
    img_geo = imread(outf1);
    figure(2); imshow(img_geo)

    % Lectura de las homografias de las imagenes
    filename3=('C:\Users\payo\Dropbox\COPY\siftDemoV4\Teste_paper_16_02_2018\correcciones_IEEE_ACCESS_28_09_2019\img_dif_sensor\homografias');
    filename33=sprintf('matriz_projective_%d.mat',ii);
    outf3=fullfile(filename3,filename33);
    hmatrixtemp= load(outf3);
    hmatrixtemp= hmatrixtemp.matriz_projective';    
    tform = affine2d(hmatrixtemp');
    
    dim  = size(img_geo);
    img_vantz = imwarp(img_vant, tform);
    figure(3); imshow(img_vantz);

    img_geo_copy = img_geo;
    img_vant_copy = img_vant;
     
    figure(4);  imshowpair(img_geo_copy,img_vant_copy,'montage','Scaling',...
    'joint');
    
    img_vant_match = imhistmatch(img_vant_copy,img_geo_copy,256); 

    if size(img_vant_match,3)==3
        img_vant_match=rgb2gray(img_vant_match);
    else
        img_vant_match=img_vant_copy
    end

    img_vant_match_eq = histeq(img_vant_match);

    
    if size(img_geo_copy,3)==3
        img_geo_gray=rgb2gray(img_geo_copy);
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

    img_geo_novo = medfilt2(img_geo_eq_bl ,[3 3]);
    img_vant_novo = medfilt2(img_vant_match_eq_bl_gray,[3 3]);
    figure(5);  imshowpair(img_geo_novo,img_vant_novo,'montage','Scaling','joint');
    pause
end
close all;
clc