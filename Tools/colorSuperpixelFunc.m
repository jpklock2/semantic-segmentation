function colorSuperpixelFunc(rgbImage, classes, idx, parameters)

N = length(classes);
outputImage = zeros(size(rgbImage),'like',double(rgbImage));
numRows = size(outputImage,1);
numCols = size(outputImage,2);
for labelVal = 1:N
    redIdx = idx{labelVal};
    greenIdx = idx{labelVal}+numRows*numCols;
    blueIdx = idx{labelVal}+2*numRows*numCols;
    outputImage(redIdx) = parameters.meanRed(classes(labelVal));
    outputImage(greenIdx) = parameters.meanGreen(classes(labelVal));
    outputImage(blueIdx) = parameters.meanBlue(classes(labelVal));
end
figure; montage({rgbImage,outputImage});

end
