%%%
% Por: Carlos Anderson (calicinio@gmail.com)
% 
% in:
%  image1 - Imagem que será o target da correlação
%  image2 - Imagem que será o template da correlação
%
% out:
%  corrMat - Matriz de correlação entre as imagens
%  centro - Coordenada do centro da image1 com maior correlação 
%
%%%

%Realiza o calculo da matriz de correlação entre duas imagens.
function [corrMat, centro]=correlacao(image1,image2)

% Converte as imagens para tons de cinza, 
% caso a mesma ainda não tenha sido convertida
if size(image1,3)==3
    image1=rgb2gray(image1);
end
if size(image2,3)==3
    image2=rgb2gray(image2);
end

% Filtro a ser utilizado.
%filtro = 'canny';

% Aplicação do filtro
%image1_bordas = edge(image1, filtro);
%image2_bordas = edge(image2, filtro);

if size(image1)>size(image2_bordas)
    Target=image1_bordas;
    Template=image2_bordas;
else
    Target=image2_bordas;
    Template=image1_bordas;
end

% Extrai as medidas das imagens
[r1,c1]=size(Target);
[r2,c2]=size(Template);

% figure
% imshow(Target)
% 
% figure
% imshow(Template)

% Calcula e monta matriz de correlação
image22=Template-mean(mean(Template));
corrMat=[];
for i=1:(r1-r2+1)
    for j=1:(c1-c2+1)
        Nimage=Target(i:i+r2-1,j:j+c2-1);
        Nimage=Nimage-mean(mean(Nimage));
        corr=sum(sum(Nimage.*image22));
        corrMat(i,j)=corr;
    end 
end


% Extrai as coordenadas de centro o ponto 
% de maior correlação entre as imagens
[r,c]=max(corrMat);
[r3,c3]=max(max(corrMat));

i=c(c3);
j=c3;

ptCenterX = floor(i + (r2/2));
ptCenterY = floor(j + (c2/2));

centro = [ptCenterX ptCenterY];
