if printResults
fprintf('\nSegmenting superpixels...\n');
tic;
end
 [L,N] = superpixels(rgbImage, K, 'Method', 'slic', 'Compactness', 20, 'NumIterations', 10);
%  [L,N] = superpixels(rgbImage, K, 'Method', 'slic0');
 idx = label2idx(L);
%      fprintf('\nTempo para segmentaaao de K = %f superpixels: %f s\n', K, toc);
a = 0;
if a
    BW = boundarymask(L);
    figure; montage({rgbImage, imoverlay(rgbImage,BW,'cyan')});
end
%     subplot(2,4,3);
%     if ishsv
%         rgbPlotImage = hsv2rgb(rgbImage);
%         imshow(imoverlay(rgbPlotImage,BW,'cyan'),'InitialMagnification',67);
%     else
%         imshow(imoverlay(rgbImage,BW,'cyan'),'InitialMagnification',67);
%     end
%     title('Super Pixel Segments')
%  end
if printResults
fprintf('Execution time superpixel segmentation: %f s\n', toc);
end