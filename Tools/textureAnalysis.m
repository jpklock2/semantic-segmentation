function textureAnalysis(parameters, centroids, L, idx, classes, rgbImage)

%% Load colors and init labels
loadColors;
Nc = length(unique(classes));
labelsArr = cell(1, Nc);
for i = 1:length(labelsArr)
  labelsArr{i} = ['c' num2str(i)];    
end
myColors(1,:) = [];

%% Get Membership Functions

inFis = parameters.fis;
cent = zeros(size(inFis.Inputs, 2), size(inFis.Inputs(1).MembershipFunctions, 2));
sig = zeros(size(inFis.Inputs, 2), size(inFis.Inputs(1).MembershipFunctions, 2));
x = -0.5:0.0001:0.5;
mfsInp = {[]};
for j=1:size(inFis.Inputs(1).MembershipFunctions, 2)
    figure(2); subplot(2,4,j); hold on;
    mfs = [];
	for i=1:size(inFis.Inputs, 2)
        if i == 1 || i == 2; plot(x, gaussmf(x, inFis.Inputs(i).MembershipFunctions(j).Parameters)); end
        mfs = [mfs; gaussmf(x, inFis.Inputs(i).MembershipFunctions(j).Parameters)];
        cent(i,j) = inFis.Inputs(i).MembershipFunctions(j).Parameters(2);
        sig(i,j) = inFis.Inputs(i).MembershipFunctions(j).Parameters(1);
	end
    title(['Centroid ' num2str(j)], 'Color', myColors(j, :));
    figure(3); subplot(2,4,j); plot(x, max(mfs(1:2, :))); title(['Centroid ' num2str(j)]);
    mfsInp = [mfsInp; mfs];
end

figure(2);
myLegs = {[]};
for i = 1:Nc
    myLegs = [myLegs; ['MF' num2str(i)]];
end
myLegs = myLegs(2:end);
myLeg = legend(myLegs, 'Orientation','horizontal');
newPosition = [0.4 0.01 0.2 0.05];
newUnits = 'normalized';
set(myLeg,'Position', newPosition,'Units', newUnits);

mfsInp = mfsInp(2:end);
porcs = [0.3, 0.7];
figure;
for i = 1:2
mf = max(porcs(i)*max(mfsInp{3}),(1-porcs(i))*max(mfsInp{7}));
subplot(1,2,i); plot(x, mf); hold on; title(['Ex' num2str(i) ': ' num2str(100*porcs(i)) '% c3 - ' num2str(100*(1-porcs(i))) '% c7']);

xCentroid1 = defuzz(x,mf,'centroid');
idxC = find(x >= xCentroid1);
kValue = mf(idxC(1));

xCentroid2 = defuzz(x,mf,'bisector');
idxC2 = find(x >= xCentroid2);
kValue2 = mf(idxC2(1));

hCentroid = line([xCentroid1 xCentroid1],[min(mf) max(mf)],'Color','k');
hCentroid = line([xCentroid2 xCentroid2],[min(mf) max(mf)],'Color','r');
end
myLeg = legend({'MF', 'Defuzzy Centroid', 'Defuzzy Bisector'}, 'Orientation','horizontal');
newPosition = [0.4 0.01 0.2 0.05];
newUnits = 'normalized';
set(myLeg,'Position', newPosition,'Units', newUnits);

%% Get centroids
centroidsReconstructed = ((parameters.pcaCoeffs(:,1:parameters.pcaN)*centroids')+parameters.pcaMean')';
textureCentroids = centroidsReconstructed(:, 85:end);

%% Box Plots

figure;
for i = 1:Nc
    subplot(2,ceil(Nc/2),i); boxplot(textureCentroids(i, :)); title(['Class ' num2str(i)], 'Color', myColors(i, :));
end
imlegend(myColors(1:Nc, :), labelsArr, 0.01, 0.05);

%% Plot texture coefficients
figure;
for i = 1:Nc
    subplot(2,ceil(Nc/2),i); plot(textureCentroids(i, :)); title(['Class ' num2str(i)], 'Color', myColors(i, :));
    hold on; plot(1:size(textureCentroids, 2), ones(1, size(textureCentroids, 2))*mean(textureCentroids(i, :)), 'r');
end
myLeg = legend({'Coefficients', 'Mean'}, 'Orientation','horizontal');
newPosition = [0.4 0.01 0.2 0.05];
newUnits = 'normalized';
set(myLeg,'Position', newPosition,'Units', newUnits);
suptitle('Centroids Texture Coefficients');

%% Plot data visualization (PCA)
[coeff,~,~,~,explained,mu] = pca(textureCentroids);
se = explained(1) + explained(2);
pcaC = (coeff(:,1:2)'*(textureCentroids-mu)')';
figure; subplot(121); hold on;
for i = 1:size(pcaC, 1)
    plot(pcaC(i, 1), pcaC(i, 2), '*', 'Color', myColors(i, :));
end
title(['First 2 PCA components (' num2str(round(se, 2)) '% variance)']);
axis auto;

%% Plot data visualization (t-SNE)
tC = tsne(textureCentroids);
subplot(122); hold on;
for i = 1:size(tC, 1)
    plot(tC(i, 1), tC(i, 2), '*', 'Color', myColors(i, :));
end
title('t-SNE plot');
axis auto;
suptitle('Components Comparison');

myLeg = legend(labelsArr, 'Orientation','horizontal');
newPosition = [0.4 0.01 0.2 0.05];
newUnits = 'normalized';
set(myLeg,'Position', newPosition,'Units', newUnits);

%% Plot distance matrix
distM = pdist2(textureCentroids, textureCentroids, 'correlation');
figure; imagesc(distM);
title('Correlation Between Centroids')
M=size(distM);
M2=distM;
mat_string=reshape(arrayfun(@(x) strtrim(cellstr(sprintf('%2.2f',x))),M2),M);
X = 1:8;
Y = 1:8;
for i = 1:M(1)
    for j = 1:M(2)
        if (distM(i,j) < 1)
            text(X(j),Y(i),mat_string{i,j},'horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','white');
        else
            text(X(j),Y(i),mat_string{i,j},'horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','black');
        end
    end
end

%% Get class image
classesIdx = cell(Nc, 1);
numRows = size(L,1);
numCols = size(L,2);
colorLmask = zeros(size(L, 1), size(L, 2), 3);
for i = 1:length(idx)
    colorLmask(idx{i}) = myColors(classes(i), 1);
    colorLmask(idx{i}+numRows*numCols) = myColors(classes(i), 2);
    colorLmask(idx{i}+2*numRows*numCols) = myColors(classes(i), 3);
end

%% Plot class image montage
if ~exist('rgbImage','var')
    figure; imshow(colorLmask);
    imlegend(myColors(1:Nc, :), labelsArr, 0.02);
    suptitle('Image Classes');
else
    outputImage = zeros(size(rgbImage),'like',double(rgbImage));
    for labelVal = 1:length(classes)
        outputImage(idx{labelVal}) = parameters.meanRed(classes(labelVal));
        outputImage(idx{labelVal}+numRows*numCols) = parameters.meanGreen(classes(labelVal));
        outputImage(idx{labelVal}+2*numRows*numCols) = parameters.meanBlue(classes(labelVal));
    end
    figure; montage({rgbImage,outputImage,colorLmask}, 'Size', [1 3]);
    imlegend(myColors(1:Nc, :), labelsArr, 0.02, 0.1);
    suptitle('Image Classes');
end

%% Plot final result text 1
sigC = 1./(1+exp(-mean(textureCentroids, 2)));
hF=figure; set(hF,'color',[1 1 1])
hA=axes; set(hA,'color',[1 1 1],'visible','off')
labelsArrSig = cell(Nc, 1);
text(0.5,0.9,'Sigmoid( mean(centroids) )','horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','black');
for i = 1:length(labelsArrSig)
    text(0.1+(0.1*(i-1)),0.7,num2str(round(sigC(i), 3)),'color',myColors(i, :));
end
text(0.5,0.5,'Mapping to 0.5 - 1.0','horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','black');
sigC = (sigC+1)/2;
for i = 1:length(labelsArrSig)
    text(0.1+(0.1*(i-1)),0.3,num2str(round(sigC(i), 3)),'color',myColors(i, :));
end
imlegend(myColors(1:Nc, :), labelsArr, 0.1, 0.1);

%% Plot final result text 2
sigC = 1./(1+exp(-mean(textureCentroids, 2)./(1*std(textureCentroids, [], 2))));
hF=figure; set(hF,'color',[1 1 1])
hA=axes; set(hA,'color',[1 1 1],'visible','off')
labelsArrSig = cell(Nc, 1);
text(0.5,0.9,'Sigmoid( mean(centroids) / std(centroids) )','horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','black');
for i = 1:length(labelsArrSig)
    text(0.1+(0.1*(i-1)),0.7,num2str(round(sigC(i), 3)),'color',myColors(i, :));
end
text(0.5,0.5,'Mapping to 0.5 - 1.0','horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','black');
sigC = (sigC+1)/2;
for i = 1:length(labelsArrSig)
    text(0.1+(0.1*(i-1)),0.3,num2str(round(sigC(i), 3)),'color',myColors(i, :));
end
imlegend(myColors(1:Nc, :), labelsArr, 0.1, 1.1);

%% Plot final result text 3
sigC = 1./(1+exp(-mean(textureCentroids, 2)./(3*std(textureCentroids, [], 2))));
hF=figure; set(hF,'color',[1 1 1])
hA=axes; set(hA,'color',[1 1 1],'visible','off')
labelsArrSig = cell(Nc, 1);
text(0.5,0.9,'Sigmoid( mean(centroids) / 3*std(centroids) )','horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','black');
for i = 1:length(labelsArrSig)
    text(0.1+(0.1*(i-1)),0.7,num2str(round(sigC(i), 3)),'color',myColors(i, :));
end
text(0.5,0.5,'Mapping to 0.5 - 1.0','horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','black');
sigC = (sigC+1)/2;
for i = 1:length(labelsArrSig)
    text(0.1+(0.1*(i-1)),0.3,num2str(round(sigC(i), 3)),'color',myColors(i, :));
end
imlegend(myColors(1:Nc, :), labelsArr, 0.1, 0.1);

%% Plot final result text 4
sigC = 1./(1+exp(-mean(textureCentroids, 2)./(3*mean(std(textureCentroids, [], 2)))));
hF=figure; set(hF,'color',[1 1 1])
hA=axes; set(hA,'color',[1 1 1],'visible','off')
labelsArrSig = cell(Nc, 1);
text(0.5,0.9,'Sigmoid( mean(centroids) / 3*mean(std(centroids)) )','horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','black');
for i = 1:length(labelsArrSig)
    text(0.1+(0.1*(i-1)),0.7,num2str(round(sigC(i), 3)),'color',myColors(i, :));
end
text(0.5,0.5,'Mapping to 0.5 - 1.0','horizontalalignment','center','verticalalignment','middle','fontsize',15,'color','black');
sigC = (sigC+1)/2;
for i = 1:length(labelsArrSig)
    text(0.1+(0.1*(i-1)),0.3,num2str(round(sigC(i), 3)),'color',myColors(i, :));
end
imlegend(myColors(1:Nc, :), labelsArr, 0.1, 0.1);

end

%% Aux function
function imlegend(colorArr, labelsArr, yp, ys)
hold on;
for ii = 1:length(labelsArr)
  hidden_h(ii) = scatter([],[],1, colorArr(ii,:), 'filled', 'DisplayName', labelsArr{ii});
end
hold off;
% uistack(hidden_h, 'bottom');
% legend(hidden_h, labelsArr);
% legend();

myLeg = legend(labelsArr, 'Orientation','horizontal');
newPosition = [0.4 yp 0.2 ys];
newUnits = 'normalized';
set(myLeg,'Position', newPosition,'Units', newUnits);
end