% seleciona pixels
close all

load('./Classes/myClasses.mat');
h4 = figure(3);
BW2 = boundarymask(L2);
h5 = imshow(imoverlay(mapImage,BW2,'cyan'),'InitialMagnification',100);
hold on;
h6 = imshow(classesImage);
h6.AlphaData = 0.35;
h7 = figure(4);
h8 = imshow(imoverlay(mapImage,BW2,'cyan'),'InitialMagnification',100);

if (useCropped)
    eval(['load(''./Classes/im' num2str(m) '.mat'')']);
    h9 = figure(5);
    BW3 = boundarymask(L2);
    h10 = imshow(imoverlay(rgbImageCrop,BW3,'cyan'),'InitialMagnification',100);
    hold on;
    h11 = imshow(tempImage);
    h11.AlphaData = 0.75;
    h12 = figure(6);
    h13 = imshow(imoverlay(rgbImageCrop,BW3,'cyan'),'InitialMagnification',100);
end

myClasses = {};
try
    if (useCropped)
        eval(['load(''./Classes/cropped_im' num2str(m) '.mat'')']);
    else
        eval(['load(''./Classes/im' num2str(m) '.mat'')']);
    end
catch
end
loaded = 1;
if isempty(myClasses)
    loaded = 0;
    % load('tempSelect.mat') % ativar manuamente se travou ou parou no meio
    tempImage = rgbImage;
end
numRows = size(rgbImage,1);
numCols = size(rgbImage,2);
counter = 0;
idx = label2idx(L);
BW = boundarymask(L);
h1 = figure(1); hold on;
h2 = imshow(imoverlay(rgbImage,BW,'cyan'),'InitialMagnification',100);
truesize
% imshow(rgbImageTemp);
L2 = L;
N2 = N;

for class = 1:length(myColors)

    if loaded && class <= length(myClasses)
        spClasses = myClasses{class};
    else
        spClasses = [];
    end
    h3 = figure(2); hold on;
    imshow(imoverlay(tempImage,BW,'cyan'),'InitialMagnification',100);
    title(['\fontsize{20} Class ' num2str(class) ' with color ' ...
           '\color[rgb]{' num2str(myColors(class, 1)) ' ' num2str(myColors(class, 2))...
           ' ' num2str(myColors(class, 3)) '}' myColorsNmes{class}]);
    truesize
    
    while ishandle(h3)
        try
            [x, y] = ginput(1);
            spNumber = L(round(y), round(x));
            spClasses = [spClasses; spNumber];
            redIdx = idx{spNumber};
            greenIdx = idx{spNumber}+numRows*numCols;
            blueIdx = idx{spNumber}+2*numRows*numCols;
            tempImage(redIdx) = myColors(class, 1);
            tempImage(greenIdx) = myColors(class, 2);
            tempImage(blueIdx) = myColors(class, 3);
            counter = counter+1;
            if mod(counter, 20) == 0
                save tempSelect.mat tempImage L N
%                 eval(['save ./Classes/im' num2str(m) '.mat tempImage L N myClasses']);
            end
            imshow(imoverlay(tempImage,BW,'cyan'),'InitialMagnification',100);
            truesize
        catch
            
        end
        
    end
    save tempSelect.mat tempImage tempImage L N
    myClasses{class} = spClasses;
    
    if (useCropped)
        eval(['save ./Classes/cropped_im' num2str(m) '.mat tempImage L2 N2 myClasses']);
    else
        eval(['save ./Classes/im' num2str(m) '.mat tempImage L2 N2 myClasses']);
    end
end
