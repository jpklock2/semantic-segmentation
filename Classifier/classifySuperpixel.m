switch myClassifier

    case 1
        % MSE
        classesP = [];
        for p = 1:N
            cla = 0; error = Inf;
            for c = 1:size(centroidsFinal, 1)
                if error > immse(pixels(p, :), centroidsFinal(c, :))
                    error = immse(pixels(p, :), centroidsFinal(c, :));
                    cla = c;
                end
            end
            classesP = [classesP; cla];
        end
        dbg = 1;
        if (myClusters == 0)
            accuracy = sum(classesP == classes)/length(classes);
    %         accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
            predAcc = predominanceAccuracy(classes, classesP);
    %         predominance = histc(classesP, unique(classesP))./length(classesP);
            myAccs = [myAccs; accuracy predAcc mean([accuracy predAcc])];
    %         myAccs = [myAccs; accuracy accuracy2 sum(predominance > 0.2) sum(predominance > 0.3) sum(predominance > 0.4) sum(predominance > 0.5)];
    %         fprintf('Accuracy = %.2f%% [%d/%d]\n', 100*accuracy, sum(classesP == classes), length(classes));
            fprintf('Accuracy = %.2f%% [%d/%d]\n\n', 100*accuracy, sum(classesP == classes), length(classes));
        end
        classes = classesP;
%         colorSuperpixel;
    case 2
        % SVM
%         addpath('E:\PJ\libsvm-3.24\matlab');
        [classesP, accuracy, prob_estimates] = svmpredict(classes, pixels, parameters.model, '-q');
        if (myClusters == 0)
    %         myAccs = [myAccs; accuracy(1)];
    %         accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
            accuracy = sum(classes == classesP ) / length(classes);
            predAcc = predominanceAccuracy(classes, classesP);
            myAccs = [myAccs; accuracy predAcc mean([accuracy predAcc])];
            fprintf('Accuracy = %.2f%% [%d/%d]\n\n', 100*accuracy, sum(classesP == classes), length(classes));
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
        YsOut = zeros(length(classes), 1);
        for k=1:length(classes)
            for o=1:nO
                [YsOut(k, o),y,w,b] = saida(pixels(k,:),parameters.p(:,:,o),parameters.q,parameters.s,parameters.c,nRe,nMf,1);
            end
            YsOut(k,:) = softmax(YsOut(k,:)')';
        end
%         YsOut = round(YsOut);
%         YsOut(YsOut == 0) = 1;
%         YsOut(YsOut == KF+1) = KF-1;
%         classes = YsOut;
        [~, idxC] = sort(YsOut, 2);
        classesP = idxC(:, end);
        if (myClusters == 0)
    %         accuracy = sum(YsOut == classes | YsOut == classes+1 | YsOut == classes-1) / length(pixels);
    %         accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
            accuracy = sum(classes == classesP ) / length(classes);
            predAcc = predominanceAccuracy(classes, classesP);
            myAccs = [myAccs; accuracy predAcc mean([accuracy predAcc])];
    %         myAccs = [myAccs; accuracy accuracy2 sum(predominance > 0.2) sum(predominance > 0.3) sum(predominance > 0.4) sum(predominance > 0.5)];
            fprintf('Accuracy = %.2f%% [%d/%d]\n\n', 100*accuracy, sum(classesP == classes), length(classes));
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
        for p = 1:N
            [~, cla] = min(pdist2(centroidsFinal, pixels(p, :), 'cosine'));
            classesP = [classesP; cla];
        end
        if (myClusters == 0)
            accuracy = sum(classesP == classes)/length(classes);
    %         accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
            predAcc = predominanceAccuracy(classes, classesP);
            myAccs = [myAccs; accuracy predAcc mean([accuracy predAcc])];
            fprintf('Accuracy = %.2f%% [%d/%d]\n', 100*accuracy, sum(classesP == classes), length(classes));
        end
        classes = classesP;

end