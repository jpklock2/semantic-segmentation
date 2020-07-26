% limpando telas
clear all;
close all;
clc;
warning ('off','all');
format long g

addpath(genpath('.'));

listing = dir('./Images/Test/Original');
imageNames = [{}];
imageNames{1} = 'Mosaicoo.tif'; % mudar para pegar automatico eventualmente
for i=1:length(listing)
    if contains(lower(listing(i).name), '.jpg')
        imageNames = [imageNames; listing(i).name];
    end
end
% imageNamesCell = cellstr(imageNames);
