
if strcmp(dataset, 'SP')
    imagesPath = 'Images/Data_BR_SP/UAV_img_partial';
    listing = dir(['./' imagesPath]);
    imageNames = [{}];
    imageNames{1} = 'Mosaicoo.tif';
    for i=1:length(listing)
        if contains(lower(listing(i).name), '.jpg')
            imageNames = [imageNames; listing(i).name];
        end
    end
    
else
    imagesPath = 'Images/Data_Suecia/UAV_img';
    listing = dir(['./' imagesPath]);
    imageNames = [{}];
    imageNames{1} = 'ortho_19comp.jpg';
    for i=1:length(listing)
        if contains(lower(listing(i).name), '.pgm')
            imageNames = [imageNames; listing(i).name];
        end
    end
end