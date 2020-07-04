fprintf('\nSegmenting superpixels...\n');
tic;
 [L,N] = superpixels(rgbImage, K, 'Method', 'slic');
%      [L,N] = superpixels(rgbImage, K, 'Method', 'slic0');
 idx = label2idx(L);
%      fprintf('\nTempo para segmentaaao de K = %f superpixels: %f s\n', K, toc);
a = 0; 
if a
    BW = boundarymask(L);
    figure; imshow(imoverlay(rgbImage,BW,'cyan'),'InitialMagnification',67);
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
fprintf('Execution time superpixel segmentation: %f s\n', toc);
