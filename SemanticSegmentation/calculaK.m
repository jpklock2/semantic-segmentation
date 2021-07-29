if printResults
fprintf('\nCalculating best K...\n');
tic;
end
%% K otimo por joelho da funaao objetivo
if (bestK == 0)
    tic;
    jotas = [];
    for k = 2:30
        jotasInd = [];
        for k2 = 1:10
            [~, ~, Jc, ~] = MyFuzzyMeans_opt(pixels, k);
            jotasInd = [jotasInd; Jc(end)];
%         [~, Utemp, J, centroids] = MyFuzzyMeans_opt(pixels, k);
%         [~, idx2] = sort(Utemp, 2);
%         jotas = [jotas; J(end)];
        end
        jotas = [jotas; median(jotasInd)];
    end
    KF = knee_pt(jotas)+2;
    if KF <= 0
        KF = 1;
    end
%     KF = knee_pt(jotas)+2;

%     if plotsIM
%         [outputImage] = evalFunction(pixels, KF, idx, A, N);
%         subplot(2,4,8);
%         if ishsv
%             outputImage = hsv2rgb(outputImage);
%         end
%         imshow(outputImage,'InitialMagnification',67);
%         title(['Objective Function K=' num2str(KF)])
%     end
end

%% K atimo por CalinskiHarabasz
if (bestK == 1)
    tic; % inicia contagem do tempo
    klist=2:30;%the number of clusters you want to try
    myfunc = @(X,K)(MyFuzzyMeans_opt(X, K));
    eva = evalclusters(pixels,myfunc,'CalinskiHarabasz','klist',klist);
    fprintf('\nTempo para K atimo = %f por CalinsHarabasz: %f s', eva.OptimalK, toc);
    [outputImage] = evalFunction(pixels, eva.OptimalK, idx, A, N);

    if plotsIM
        subplot(2,4,4);
        if ishsv
            outputImage = hsv2rgb(outputImage);
        end
        imshow(outputImage,'InitialMagnification',67);
        title(['CalinskiHarabasz K=' num2str(eva.OptimalK)])
    end
end

%% K atimo por Silhouette
if (bestK == 2)
    tic;
    klist=2:30;%the number of clusters you want to try
    myfunc = @(X,K)(MyFuzzyMeans_opt(X, K));
    eva = evalclusters(pixels,myfunc,'Silhouette','klist',klist);
    fprintf('\nTempo para K atimo = %f por Silhouette: %f s', eva.OptimalK, toc);
    [outputImage] = evalFunction(pixels, eva.OptimalK, idx, A, N);

    if plotsIM
        subplot(2,4,5);
        if ishsv
            outputImage = hsv2rgb(outputImage);
        end
        imshow(outputImage,'InitialMagnification',67);
        title(['Silhouette K=' num2str(eva.OptimalK)])
    end
end

%% K atimo por DaviesBouldin
if (bestK == 3)
    tic;
    klist=2:30;%the number of clusters you want to try
    myfunc = @(X,K)(MyFuzzyMeans_opt(X, K));
    eva = evalclusters(pixels,myfunc,'DaviesBouldin','klist',klist);
    fprintf('\nTempo para K atimo = %f por DaviesBouldin: %f s', eva.OptimalK, toc);
    [outputImage] = evalFunction(pixels, eva.OptimalK, idx, A, N);

    if plotsIM
        subplot(2,4,6);
        if ishsv
            outputImage = hsv2rgb(outputImage);
        end
        imshow(outputImage,'InitialMagnification',67);
        title(['DaviesBouldin K=' num2str(eva.OptimalK)])
    end
end

%% K atimo por GAP
if (bestK == 4)
    tic;
    klist=2:30;%the number of clusters you want to try
    myfunc = @(X,K)(MyFuzzyMeans_opt(X, K));
    eva = evalclusters(pixels,myfunc,'GAP','klist',klist);
    fprintf('\nTempo para K atimo = %f por GAP: %f s\n', eva.OptimalK, toc);
    [outputImage] = evalFunction(pixels, eva.OptimalK, idx, A, N);

    if plotsIM
        subplot(2,4,6);
        imshow(outputImage,'InitialMagnification',67);
        title(['GAP K=' num2str(eva.OptimalK)])
    end
end

%% K otimo por silhueta alternativo
if (bestK == 5)
    tic;
    sil = [];
    for k=3:20
        [~, Utemp, J, centroids] = MyFuzzyMeans_opt(pixels, k);
        [~, idx2] = sort(Utemp, 2);
        sil = [sil; mean(silhouette(pixels, idx2(:, end)))];
    end
%     K2 = knee_pt(sil)+3;
    [~, idxS] = sort(sil);
    KF = idxS(end)+2;
    fprintf('\nTempo para K atimo = %f por Silhueta: %f s', KF, toc);
    [outputImage] = evalFunction(pixels, KF, idx, A, N);

    if plotsIM
        subplot(2,4,7);
        if ishsv
            outputImage = hsv2rgb(outputImage);
        end
        imshow(outputImage,'InitialMagnification',67);
        title(['Alternative Silhouette K=' num2str(KF)])
    end
end
if printResults
fprintf('Best K = %d', KF);
fprintf('\nExecution time for calculating best K: %f s\n', toc);
end