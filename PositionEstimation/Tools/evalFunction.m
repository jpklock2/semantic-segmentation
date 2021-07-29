 function [outputImage] = evalFunction(classes, K, idx, A, N)
 
    plots = 0; % plot em 2 dimensaes
    plots3 = 0; % plot em 3 dimensaes
    prints = 0; % print dos valores dos resultados
    
    numRows = size(A,1);
    numCols = size(A,2);
 
%     [~, U, ~, centroids] = MyFuzzyMeans_opt(pixels, K);
    
    outputImage = zeros(size(A),'like',A);
%     clus = 1:K; % declara clusters
%     [~, idx2] = sort(U, 2); % ordena matriz de probabilidades por linha (elemento)
%     classesTemp = idx2(:, end);
    for c = 1:K
        idxC = find(classes == c);
        redIdx = [];
        greenIdx  = [];
        blueIdx = [];
        for j = 1:length(idxC)
            tempRed = idx{idxC(j)};
            redIdx = [redIdx; tempRed];
            tempGreen = idx{idxC(j)}+numRows*numCols;
            greenIdx = [greenIdx; tempGreen];
            tempBlue = idx{idxC(j)}+2*numRows*numCols;
            blueIdx = [blueIdx; tempBlue];
        end
        outputImage(redIdx) = mean(A(redIdx));
        outputImage(greenIdx) = mean(A(greenIdx));
        outputImage(blueIdx) = mean(A(blueIdx));
    end
 end
  