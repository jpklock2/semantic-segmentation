function [predAcc] = predominanceAccuracy(classes, classesP)

predominance = histc(classes, unique(classes))./length(classes);
predominanceP = histc(classesP, unique(classesP))./length(classesP);
if (length(predominance) > length(predominanceP)); predominanceP(length(predominanceP)+1:length(predominance)) = 0; end
if (length(predominanceP) > length(predominance)); predominance(length(predominance)+1:length(predominanceP)) = 0; end
if (size(predominanceP, 2) > size(predominanceP, 1)); predominanceP = predominanceP'; end
if (size(predominance, 2) > size(predominance, 1)); predominance = predominance'; end
predominanceP = sort(predominanceP);
predominance = sort(predominance);
tempP = zeros(length(predominance), 1);
tempP(end) = 1;
maxD = max(sum(abs(predominance-(1/length(predominance)))), sum(abs(predominance-tempP)));
predAcc = 1-sum(abs(predominance-predominanceP))/maxD;

end
