clear
close all

lambda = 1;
iter = 4;
p = 0.8;
eps = 0.0001; 

Img = im2double((imread('flower.png')));
Smoothed = ILS_LNorm(Img, lambda, p, eps, iter);
% Smoothed = ILS_LNorm_GPU(Img, lambda, p, eps, iter);

Diff = Img - Smoothed;
ImgE = Img + 3 * Diff;

figure; imshow(Smoothed)
figure; imshow(ImgE)

