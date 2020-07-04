switch myClassifier

    case 1
        % RMSE
        classesP = [];
        for p = 1:size(centroids, 1)
            cla = 0; error = Inf;
            for c = 1:size(centroidsFinal, 1)
                if error > sqrt(immse(centroids(p, :), centroidsFinal(c, :)))
                    error = sqrt(immse(centroids(p, :), centroidsFinal(c, :)));
                    cla = c;
                end
            end
            classesP = [classesP; cla];
        end
        dbg = 1;
        classesTemp2 = classesTemp;
        for i = 1:length(pixels)
            classesTemp2(i) = classesP(classesTemp2(i));
        end
        classesP = classesTemp2;
        accuracy = sum(classesP == classes)/length(classes);
%         accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
        predAcc = predominanceAccuracy(classes, classesP);
        myAccs = [myAccs; accuracy predAcc mean([accuracy predAcc])];
%         myAccs = [myAccs; accuracy accuracy2 sum(predominance > 0.2) sum(predominance > 0.3) sum(predominance > 0.4) sum(predominance > 0.5)];
%         fprintf('Accuracy = %.2f%% [%d/%d]\n', 100*accuracy, sum(classesP == classes), length(classes));
        if (myClusters == 0)
            fprintf('Accuracy = %.2f%% [%d/%d]\n', 100*accuracy, sum(classesP == classes), length(classes));
        end
        classes = classesP;
    case 2
        % SVM
%         addpath('E:\PJ\libsvm-3.24\matlab');
        [classesP, accuracy, prob_estimates] = svmpredict(rand([size(centroids, 1) 1]), pixels, model, '-q');
        classesTemp2 = classesTemp;
        for i = 1:length(pixels)
            classesTemp2(i) = classesP(classesTemp2(i));
        end
        classesP = classesTemp2;
        accuracy = sum(classesP == classes)/length(classes);
%         accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
        predAcc = predominanceAccuracy(classes, classesP);
        myAccs = [myAccs; accuracy predAcc mean([accuracy predAcc])];
        if (myClusters == 0)
            fprintf('Accuracy = %.2f%% [%d/%d]\n', 100*accuracy, sum(classesP == classes), length(classes));
        end
        classes = classesP;
%         colorSuperpixel;
        
   case 3
        % ANFIS
%         maxIter = 100;
%        inFis = genfis3(pixelsOri, classes, 'sugeno', length(unique(classes))); % gera sistema nebuloso
%        [outFis, trainError] = anfis([pixelsOri classes], inFis, maxIter); % treina ANFIS
%         YsOut = evalfis(outFis, pixels); % gera saida para dados de validacao
        nO = KF;
        nRe = nO;
        nMf=length(pixelsOri(1,:));
        YsOut = zeros(size(centroids, 1), 1);
        for k=1:size(centroids, 1)
            for o=1:nO
                [YsOut(k, o),y,w,b] = saida(centroids(k,:),pi(:,:,o),qi,sig,cent,nRe,nMf,1);
            end
            YsOut(k,:) = softmax(YsOut(k,:)')';
        end
%         YsOut = round(YsOut);
%         YsOut(YsOut == 0) = 1;
%         YsOut(YsOut == KF+1) = KF-1;
%         classes = YsOut;
        [~, idxC] = sort(YsOut, 2);
        classesP = idxC(:, end);
        
        classesTemp2 = classesTemp;
        for i = 1:length(pixels)
            classesTemp2(i) = classesP(classesTemp2(i));
        end
        classesP = classesTemp2;
%         accuracy = sum(YsOut == classes | YsOut == classes+1 | YsOut == classes-1) / length(pixels);
%         accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
        accuracy = sum(classes == classesP ) / length(classes);
        predAcc = predominanceAccuracy(classes, classesP);
        myAccs = [myAccs; accuracy predAcc mean([accuracy predAcc])];
%         myAccs = [myAccs; accuracy accuracy2 sum(predominance > 0.2) sum(predominance > 0.3) sum(predominance > 0.4) sum(predominance > 0.5)];
        if (myClusters == 0)
            fprintf('Accuracy = %.2f%% [%d/%d]\n', 100*accuracy, sum(classesP == classes), length(classes));
        end
        classes = classesP;

%         YsOut(YsOut < 0.5) = 0; % menor que 0.5 eh classe 0
%         YsOut(YsOut >= 0.5) = 1; % maior que 0.5 eh classe 1
%         accuracyV = length(find(ydv==YsOut))/length(ydv); % calcula acuracia de validacao
%        YsOut = evalfis(centroidsFinal, outFis); % gera saida para dados de validacao
       
%        sigma = sqrt(-((x-c).^2) ./ (2*log(m))); 
%        sigma = sum(sigma) / N; 
%        mfparams = [sigma, c];
%        params = computemfparams (mftype, Xout(:,i), U(j,:)', center(j,numInp+i));

    case 4
        % Cosine
        classesP = [];
        for p = 1:size(centroids, 1)
            [~, cla] = min(pdist2(centroidsFinal, centroids(p, :), 'cosine'));
            classesP = [classesP; cla];
        end
        dbg = 1;
        classesTemp2 = classesTemp;
        for i = 1:length(pixels)
            classesTemp2(i) = classesP(classesTemp2(i));
        end
        classesP = classesTemp2;
        accuracy = sum(classesP == classes)/length(classes);
%         accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
        predAcc = predominanceAccuracy(classes, classesP);
        myAccs = [myAccs; accuracy predAcc mean([accuracy predAcc])];
        if (myClusters == 0)
            fprintf('Accuracy = %.2f%% [%d/%d]\n', 100*accuracy, sum(classesP == classes), length(classes));
        end
        classes = classesP;

end