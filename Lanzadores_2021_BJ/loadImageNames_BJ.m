
if strcmp(dataset, 'SP') || strcmp(dataset, 'SP_segmentation')
    if strcmp(dataset, 'SP')
        imagesPath = 'Images/Data_BR_SP/UAV_img_partial';
    elseif strcmp(dataset, 'SP_segmentation')
        imagesPath = 'Lanzadores_2021_BJ/Images/Data_BR_SP/UAV_img_partial';
    end
    listing = dir(['./' imagesPath]);
    imageNames = [{}];
    imageNames{1} = 'Mosaicoo.tif';
    for i=1:length(listing)
        if contains(lower(listing(i).name), '.jpg')
            imageNames = [imageNames; listing(i).name];
        end
    end
    
elseif strcmp(dataset, 'SC') || strcmp(dataset, 'SC_segmentation')
    if strcmp(dataset, 'SC')
        imagesPath = 'Images/Data_Suecia/UAV_img';
    elseif strcmp(dataset, 'SC_segmentation')
        imagesPath = 'Lanzadores_2021_BJ/Images/Data_Suecia/UAV_img';
    end
    listing = dir(['./' imagesPath]);
    imageNames = [{}];
    imageNames{1} = 'ortho_19comp.jpg';
    for i=1:length(listing)
        if contains(lower(listing(i).name), '.pgm')
            imageNames = [imageNames; listing(i).name];
        end
    end
end