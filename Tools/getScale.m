clc;
addpath(genpath('.'));

geo = imread('Images\Train\Mosaicoo.tif');

filename = ('Images\Test\Original');
listing = dir(filename);
imageNames = [{}];
for i=1:length(listing)
    if contains(lower(listing(i).name), '.jpg')
        imageNames = [imageNames; listing(i).name];
    end
end

filenameCrop = ('Images\Test\Cropped');
listing = dir(filenameCrop);
imageNamesCrop = [{}];
for i=1:length(listing)
    if contains(lower(listing(i).name), '.jpg')
        imageNamesCrop = [imageNamesCrop; listing(i).name];
    end
end

load('transformData.mat');
fprintf('\nHomography Scales 1:');
fprintf('\nScale X: %.3f   Scale Y: %.3f\n', tform.T(1,1), tform.T(2,2));

load('transformData_2.mat');
fprintf('\nHomography Scales 2:');
fprintf('\nScale X: %.3f   Scale Y: %.3f\n', tform.T(1,1), tform.T(2,2));

xScale = [];
yScale = [];
xScaleUav = [];
yScaleUav = [];
fprintf('\nCropped Image Scales:\n\n');
for i = 1:length(imageNames)
    out = fullfile(filename, strtrim(imageNames{i}));
    outCrop = fullfile([filenameCrop '\cropped_image_' num2str(i+1) '.JPG']);
    uav = imread(out);
    uavCrop = imread(outCrop);
    xScale = [xScale; size(uavCrop, 1)/size(geo, 1)];
    yScale = [yScale; size(uavCrop, 2)/size(geo, 2)];
    xScaleUav = [xScaleUav; size(uavCrop, 1)/size(uav, 1)];
    yScaleUav = [yScaleUav; size(uavCrop, 2)/size(uav, 2)];
    fprintf('Image %d - Scale X: %.3f   Scale Y: %.3f\n', i, size(uavCrop, 1)/size(uav, 1), size(uavCrop, 2)/size(uav, 2));
%     disp([round(size(uav, 1)*tform.T(1,1)) round(size(uav, 2)*tform.T(2,2))]);
%     disp([size(uavCrop, 1) size(uavCrop, 2)]);
end

fprintf('\nMean X: %.3f   Mean Y: %.3f\n\n', mean(xScaleUav), mean(yScaleUav));
