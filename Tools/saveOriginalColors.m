meanRed = [];
meanGreen = [];
meanBlue = [];
tempCentroids = zeros(KF, size(pixels, 2));
numRows = size(rgbImage,1);
numCols = size(rgbImage,2);
if isHSV
    rgbImage = hsv2rgb(rgbImage);
end
if isLab
    rgbImage = lab2rgb(rgbImage,'Out','double');
end
for c = 1:KF
     idxC = find(classes == c);
     if isempty(idxC)
        meanRed = [meanRed; 0];
        meanGreen = [meanGreen; 0];
        meanBlue = [meanBlue; 0];
        continue; 
     end
     redIdx = [];
     greenIdx  = [];
     blueIdx= [];
     for j = 1:length(idxC)
        tempRed = idx{idxC(j)};
        redIdx = [redIdx; tempRed];
        tempGreen = idx{idxC(j)}+numRows*numCols;
        greenIdx = [greenIdx; tempGreen];
        tempBlue = idx{idxC(j)}+2*numRows*numCols;
        blueIdx = [blueIdx; tempBlue];
     end
     tempCentroids(c, :) = mean(pixelsOri(idxC, :));
     meanRed = [meanRed; mean(rgbImage(redIdx))];
     meanGreen = [meanGreen; mean(rgbImage(greenIdx))];
     meanBlue = [meanBlue; mean(rgbImage(blueIdx))];
end
