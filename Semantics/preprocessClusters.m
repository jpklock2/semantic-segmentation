fprintf('\nDefining Cluster Number for Validation Image...\n');
tic;
jotas = [];
for k=2:30
    [~, Uc, Jc, centroids] = MyFuzzyMeans_opt(pixels, k);
    jotas = [jotas; Jc(end)];
end
Kc = knee_pt(jotas)-2;
% Kc = 7;
% for Kc = 20:-1:10
fprintf('Cluster Number = %d\n', Kc);
[~, Uc, ~, centroids] = MyFuzzyMeans_opt(pixels, Kc);
[~, idxUc] = sort(Uc, 2);
classesTemp = idxUc(:, end);
fprintf('Execution time for defining clusters: %f s\n', toc);
if plotsCompare
    [outputImage] = evalFunction(classesTemp, Kc, idx, rgbImage, N);
    figure; montage({rgbImage,outputImage});
end
% end
% classes = classesTemp;
% colorSuperpixel;
dbg = 1;