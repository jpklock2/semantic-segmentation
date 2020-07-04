% separa pixels
if isHSV
    rgbImage = hsv2rgb(rgbImage);
end
if isLab
    rgbImage = lab2rgb(rgbImage,'Out','double');
end
outputImage = zeros(size(rgbImage),'like',rgbImage);
numRows = size(outputImage,1);
numCols = size(outputImage,2);
for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    outputImage(redIdx) = meanRed(classes(labelVal));
    outputImage(greenIdx) = meanGreen(classes(labelVal));
    outputImage(blueIdx) = meanBlue(classes(labelVal));
end
% figure; imshow(outputImage);
figure; montage({rgbImage,outputImage})
