imageNames = [{}];
if strcmp(dataset, 'BR')
    imagesPath = 'Images/Dataset_BR/UAV_img';
    imageNames{1} = 'map.tif';
elseif strcmp(dataset, 'SW')
    imagesPath = 'Images/Dataset_SW/UAV_img';
    imageNames{1} = 'map.jpg';
end
listing = dir(['./' imagesPath]);
for i = 3:length(listing)
    imageNames = [imageNames; listing(i).name];
end