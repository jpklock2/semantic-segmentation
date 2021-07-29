accCombined = [];

%     case 1
% MSE
classesP = [];
for p = 1:size(centroids, 1)
    cla = 0; error = Inf;
    for c = 1:size(centroidsFinal, 1)
        if error > immse(centroids(p, :), centroidsFinal(c, :))
            error = immse(centroids(p, :), centroidsFinal(c, :));
            cla = c;
        end
    end
    classesP = [classesP; cla];
end
classesTemp2 = classesTemp;
for i = 1:length(pixels)
    classesTemp2(i) = classesP(classesTemp2(i));
end
classesP = classesTemp2;
accuracy = sum(classesP == classes)/length(classes);
% accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
predAcc = predominanceAccuracy(classes, classesP);
accCombined = [accCombined accuracy predAcc mean([accuracy predAcc])];

%     case 2
% SVM
[classesP, accuracy, prob_estimates] = svmpredict(rand([size(centroids, 1) 1]), centroids, parameters.model, '-q');
classesTemp2 = classesTemp;
for i = 1:length(pixels)
    classesTemp2(i) = classesP(classesTemp2(i));
end
classesP = classesTemp2;
accuracy = sum(classesP == classes)/length(classes);
% accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
predAcc = predominanceAccuracy(classes, classesP);
accCombined = [accCombined accuracy predAcc mean([accuracy predAcc])];

%    case 3
% ANFIS
nO = KF;
nRe = nO;
nMf=length(pixelsOri(1,:));
YsOut = zeros(size(centroids, 1), 1);
for k=1:size(centroids, 1)
    for o=1:nO
        [YsOut(k, o),y,w,b] = saida(centroids(k,:),parameters.p(:,:,o),parameters.q,parameters.s,parameters.c,nRe,nMf,1);
    end
    YsOut(k,:) = softmax(YsOut(k,:)')';
end
[~, idxC] = sort(YsOut, 2);
classesP = idxC(:, end);
classesTemp2 = classesTemp;
for i = 1:length(pixels)
    classesTemp2(i) = classesP(classesTemp2(i));
end
classesP = classesTemp2;
% accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
accuracy = sum(classes == classesP ) / length(classes);
predAcc = predominanceAccuracy(classes, classesP);
accCombined = [accCombined accuracy predAcc mean([accuracy predAcc])];

%      case 4
% Cosine
classesP = [];
for p = 1:size(centroids, 1)
    [~, cla] = min(pdist2(centroidsFinal, centroids(p, :), 'cosine'));
    classesP = [classesP; cla];
end
classesTemp2 = classesTemp;
for i = 1:length(pixels)
    classesTemp2(i) = classesP(classesTemp2(i));
end
classesP = classesTemp2;
accuracy = sum(classesP == classes)/length(classes);
% accuracy2 = sum(classes == classesP | classes == classesP+1 | classes == classesP-1) / length(classes);
predAcc = predominanceAccuracy(classes, classesP);
accCombined = [accCombined accuracy predAcc mean([accuracy predAcc])];

accCombinedAll = [accCombinedAll; accCombined];
