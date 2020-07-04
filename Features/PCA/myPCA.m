if usePCA == 1 || usePCA == 2
    
fprintf('\nApplying PCA...\n');
tic;
if (m == 1)
    pixelBkpOri = pixels;
    [coeff,~,~,~,explained,mu] = pca(pixels);
    sum_explained = 0;
    idxPCA = 0;
    myVariance = 99;
    while sum_explained < myVariance
        idxPCA = idxPCA + 1;
        sum_explained = sum_explained + explained(idxPCA);
    end
    fprintf('Number of PCA componentes for %d features and %d%% of variance: %d\n', size(pixelBkpOri,2), myVariance, idxPCA);
    pixels = (coeff(:,1:idxPCA)'*(pixels-mu)')';
%             pixelsReconstructed = ((coeff(:,1:idx95)*pixelsTemp')+mu')';
else
    pixelBkp = pixels;
    if usePCA == 2
        [coeff,~,~,~,~,mu] = pca(pixels);
    end
    pixels = (coeff(:,1:idxPCA)'*(pixels-mu)')';
end
fprintf('Execution time for applying PCA: %f s\n', toc);

end