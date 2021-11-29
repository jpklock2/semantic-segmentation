function [geo_img, X1,cmap,R,bbox,cropSize,geo_mask]=get_geo_subimg(...
    lat_req, lon_req, tamx, tamy, X, cmap, R, bbox, hab, mask)

% geo_img       - imagem do banco de dados recortada com base na lat,lon.
% cmap          - mapa de cores
% R             - matriz de referencia
% BBOX          - borda
% hab           - habilitador: 0 quando é (lat,lon). 
%                 1 quando sao coordenadas espaciales (px).   

% lat_r = row ==> x no matlab (de izquierda a derecha)
% lon_r = col ==> y no matlab (de arriba para abajo)
% latlon2pix e map2pix fornecem os mesmos pixels

if hab==0
    [lat_r,lon_r]=latlon2pix(R,lat_req,lon_req);
end
if hab==1
    [lat_r,lon_r]=map2pix(R,lon_req,lat_req);
end

pos= [lon_r lat_r 30];
color={'red'};
X1= insertShape(X,'circle', pos, 'Color', color, 'LineWidth',25, 'Opacity',0.7);

%adiciona um erro para simular o sensor inercial
% rand_x = (randi([-round(tamx/4) round(tamx/4)],1,1));
% rand_y = (randi([-round(tamy/4) round(tamy/4)],1,1));
rand_x = 0;
rand_y = 0;

lat_r=round(lat_r) + rand_x;
lon_r=round(lon_r) + rand_y;

tam_img = size(X);

%%%% PARA CASO A IMAGEM SEJA RECORTADA FORA DO TAMANHO %%%%
if (lat_r - tamx) < 0 
    inicio_x = 1;
    final_x = (tamx*2)+1;
else
    inicio_x = abs(lat_r - tamx);
    final_x = abs(lat_r + tamx);
end

if (lon_r - tamy) < 0
    inicio_y = 1;
    final_y = (tamy*2)+1;    
else
    inicio_y = abs(lon_r - tamy);
    final_y = abs(lon_r + tamy);
end

if ( final_x > tam_img(1) )
    inicio_x = tam_img(1) - (tamx*2);
    final_x = tam_img(1);
end

if ( final_y > tam_img(2) )
    inicio_y = tam_img(2) - (tamy*2);
    final_y = tam_img(2);
end

res_mask = imresize(mask, [size(X, 1) size(X, 2)], 'nearest');
geo_mask = res_mask(inicio_x:final_x, inicio_y:final_y,:);

geo_img=X(inicio_x:final_x, inicio_y:final_y,:);
X1=X1(inicio_x:final_x, inicio_y:final_y,:);
[R(3,2),R(3,1)]=pix2latlon(R,inicio_x,inicio_y);
cropSize = [inicio_x inicio_y final_x final_y];
% [R(3,2),R(3,1)]=pix2map(R,inicio_x,inicio_y);
end

