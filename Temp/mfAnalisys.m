% close all
inFis = parameters.fis;
cent = zeros(size(inFis.Inputs, 2), size(inFis.Inputs(1).MembershipFunctions, 2));
sig = zeros(size(inFis.Inputs, 2), size(inFis.Inputs(1).MembershipFunctions, 2));
x = -0.5:0.0001:0.5;
mfsInp = {[]};
for j=1:size(inFis.Inputs(1).MembershipFunctions, 2)
    figure(6); subplot(2,4,j); hold on;
    mfs = [];
	for i=1:size(inFis.Inputs, 2)
        plot(x, gaussmf(x, inFis.Inputs(i).MembershipFunctions(j).Parameters));
        mfs = [mfs; gaussmf(x, inFis.Inputs(i).MembershipFunctions(j).Parameters)];
        cent(i,j) = inFis.Inputs(i).MembershipFunctions(j).Parameters(2);
        sig(i,j) = inFis.Inputs(i).MembershipFunctions(j).Parameters(1);
	end
    title(['Centroid ' num2str(j)]);
    figure(7); subplot(2,4,j); plot(x, max(mfs)); title(['Centroid ' num2str(j)]);
    mfsInp = [mfsInp; mfs];
end

mfsInp = mfsInp(2:end);
mf = max(0.7*max(mfsInp{8}),0.3*max(mfsInp{5}));
% figure; plot(max(mfsInp{1}));
% figure; plot(max(mfsInp{5}));
figure; plot(mf);

xCentroid = defuzz(x,mf,'centroid');
idx = find(x >= xCentroid);
kValue = mf(idx(1));

xCentroid = defuzz(x,mf,'bisector');
idx = find(x >= xCentroid);
kValue2 = mf(idx(1));

xCentroid = defuzz(x,mf,'mom');
idx = find(x >= xCentroid);
kValue3 = mf(idx(1));

xCentroid = defuzz(x,mf,'lom');
idx = find(x >= xCentroid);
kValue4 = mf(idx(1));

xCentroid = defuzz(x,mf,'som');
idx = find(x >= xCentroid);
kValue5 = mf(idx(1));
% Indicate the centroid defuzzification result on the original plot.
figure('Tag','defuzz')
plot(x,mf,'LineWidth',3)
h_gca = gca;
h_gca.YTick = [0 .5 1] ;
ylim([-1 1])
hCentroid = line([xCentroid xCentroid],[-0.2 1.2],'Color','k'); 
tCentroid = text(xCentroid,-0.2,' centroid','FontWeight','bold');
% test = (parameters.pcaCoeffs(85:end,1:parameters.pcaN)'*(myTerrain-parameters.pcaMean(85:end))')';
% distM = pdist2(test, centroidsFinal, 'correlation');

% myTerrain = 0.7*textureCentroids(1,:)+0.3*textureCentroids(5,:);
% test = (parameters.pcaCoeffs(85:end,1:parameters.pcaN)'*(myTerrain-parameters.pcaMean(85:end))')';
% distM = pdist2(test, centroidsFinal, 'correlation');