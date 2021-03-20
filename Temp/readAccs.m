allAccsCent = zeros(11, 12, 10);
for i = 1:10
    load(['acc_centroids_' num2str(i) '.mat']);
    allAccsCent(1:10, :, i) = accCombinedAll;
end
allAccsCentFinal = mean(allAccsCent, 3);
allAccsCentFinal(11, :) = mean(allAccsCentFinal(1:10, :));

writematrix(round(allAccsCentFinal, 2),'centroidsAcc.csv','delimiter',';');

allAccsSp = zeros(11, 12, 10);
for i = 1:10
    load(['acc_superpixels_' num2str(i) '.mat']);
    allAccsSp(1:10, :, i) = accCombinedAll;
end
allAccsSpFinal = mean(allAccsSp, 3);
allAccsSpFinal(11, :) = mean(allAccsSpFinal(1:10, :));

writematrix(round(allAccsSpFinal, 2),'superpixelsAcc.csv','delimiter',';');

dbg = 1;
