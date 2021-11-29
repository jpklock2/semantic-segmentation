if usePCA == 1 || usePCA == 2

if printResults
fprintf('\nApplying PCA...\n');
tic;
end
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
    if printResults
    fprintf('Number of PCA componentes for %d features and %d%% of variance: %d\n', size(pixelBkpOri,2), myVariance, idxPCA);
    end
    pixels = (coeff(:,1:idxPCA)'*(pixels-mu)')';
    parameters.pcaMean = mu;
    parameters.pcaCoeffs = coeff;
    parameters.pcaN = idxPCA;
    parameters.pcaExp = explained;
%     pixelsReconstructed = ((parameters.pcaCoeffs(:,1:parameters.pcaN)*pixels')+parameters.pcaMean')';
else
    pixelBkp = pixels;
    if usePCA == 2
        [coeff,~,~,~,~,mu] = pca(pixels);
        pixels = (coeff(:,1:idxPCA)'*(pixels-mu)')';
    else
        pixels = (parameters.pcaCoeffs(:,1:parameters.pcaN)'*(pixels-parameters.pcaMean)')';
    end
end
clear pixelBkpOri coeff explained mu
if printResults
fprintf('Execution time for applying PCA: %f s\n', toc);
end

end