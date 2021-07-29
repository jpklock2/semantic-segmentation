if useAll
    listing = dir('./Images/Test/Original');
    imageNames = [{}];
    imageNames{1} = 'Mosaicoo.tif'; % mudar para pegar automatico eventualmente
    for i=1:length(listing)
        if contains(lower(listing(i).name), '.jpg')
            imageNames = [imageNames; listing(i).name];
        end
    end
else
    imagesPath = 'Images/Test/Original_Dev';
%     imagesPath = 'Images/Test/Cropped/SemEqualizacao';
%     imagesPath = 'Images/Test/Cropped/ComEqualizacao';
% imagesPath = 'Imagens_artigo';
useCropped = 0;
    listing = dir(['./' imagesPath]);
    imageNames = [{}];
    imageNames{1} = 'Mosaicoo.tif'; % mudar para pegar automatico eventualmente
    for i=1:length(listing)
        if contains(lower(listing(i).name), '.jpg')
            imageNames = [imageNames; listing(i).name];
        end
    end
end
