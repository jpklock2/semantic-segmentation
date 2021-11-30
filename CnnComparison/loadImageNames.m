if strcmp(dataset, 'BR')
    imagesPath = 'PositionEstimation/Images/Dataset_BR/UAV_img';
    listing = dir(['./' imagesPath]);
    imageNames = [{}];
    imageNames{1} = 'map.tif';
    for i=1:length(listing)
        if contains(lower(listing(i).name), '.jpg')
            imageNames = [imageNames; listing(i).name];
        end
    end
    
elseif strcmp(dataset, 'SW')
    imagesPath = 'PositionEstimation/Images/Dataset_SW/UAV_img';
    listing = dir(['./' imagesPath]);
    imageNames = [{}];
    imageNames{1} = 'map.jpg';
    for i=1:length(listing)
        if contains(lower(listing(i).name), '.pgm')
            imageNames = [imageNames; listing(i).name];
        end
    end
end