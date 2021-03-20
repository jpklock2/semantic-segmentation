clc;
% load test_canny.mat;

% lambda = 5;
% iter = 5;
% p = 0.2;
% eps = 0.0001;
% A = ILS_LNorm(Template_filt, lambda, p, eps, iter);
% figure; imshow(A);

% for expe = 0:3
%     BW_Template = edge(Template_filt, 'canny', 0.7, sqrt(2));
%     BW_Template1 = Canny_v0(Template_filt, 1, 0.7, expe);
%     BW_Template2 = Canny_v0(Template_filt, 2, 0.7, expe);
%     BW_Template3 = Canny_v0(Template_filt, 3, 0.7, expe);
%     BW_Template4 = Canny_v0(Template_filt, 4, 0.7, expe);
%     BW_Template4 = Canny_v0(Template_filt, 5, 0.7, expe);
%     figure; montage({BW_Template, BW_Template4, BW_Template2, BW_Template3});
% end

% for expe = 0:3
%     BW_Target = edge(Target_filt, 'canny', 0.7*threshOut2, sqrt(2));
%     BW_Target1 = Canny_v0(Target_filt, 1, threshOut2, expe);
%     BW_Target2 = Canny_v0(Target_filt, 2, threshOut2, expe);
%     BW_Target3 = Canny_v0(Target_filt, 3, threshOut2, expe);
%     figure; subplot(221); imshow(BW_Target); subplot(222); imshow(BW_Target1);
%             subplot(223); imshow(BW_Target2); subplot(224); imshow(BW_Target3);
% end

% fprintf('\nExperiment 1: (default conv, matlab thresh)\n');
% load('experiment_04-46-47_05-27-18.mat');
% dists = [];
% for i = 1:length(log_coords1_all)
%     dists = [dists; log_coords1_all{i}(:, 5)'];
% end
% fprintf('\nMean: %.2f %.2f %2.f %.2f\n', mean(dists));
% fprintf('Median: %.2f %.2f %2.f %.2f\n', median(dists));
% fprintf('Standard deviation: %.2f %.2f %2.f %.2f\n', std(dists));
% fprintf('Max: %.2f %.2f %2.f %.2f\n', max(dists));
% fprintf('Min: %.2f %.2f %2.f %.2f\n\n', min(dists));
% 
% fprintf('\nExperiment 2: (default conv, default thresh)\n');
% load('experiment_05-33-49_06-20-13.mat');
% dists = [];
% for i = 1:length(log_coords1_all)
%     dists = [dists; log_coords1_all{i}(:, 5)'];
% end
% fprintf('\nMeans: %.2f %.2f %2.f %.2f\n', mean(dists));
% fprintf('Median: %.2f %.2f %2.f %.2f\n', median(dists));
% fprintf('Standard deviations: %.2f %.2f %2.f %.2f\n', std(dists));
% fprintf('Max: %.2f %.2f %2.f %.2f\n', max(dists));
% fprintf('Min: %.2f %.2f %2.f %.2f\n\n', min(dists));
% 
% fprintf('\nExperiment 3: (matlab conv, matlab thresh)\n');
% load('experiment_05-33-49_07-05-18.mat');
% dists = [];
% for i = 1:length(log_coords1_all)
%     dists = [dists; log_coords1_all{i}(:, 5)'];
% end
% fprintf('\nMeans: %.2f %.2f %2.f %.2f\n', mean(dists));
% fprintf('Median: %.2f %.2f %2.f %.2f\n', median(dists));
% fprintf('Standard deviations: %.2f %.2f %2.f %.2f\n', std(dists));
% fprintf('Max: %.2f %.2f %2.f %.2f\n', max(dists));
% fprintf('Min: %.2f %.2f %2.f %.2f\n\n', min(dists));
% 
% fprintf('\nExperiment 4: (matlab conv, default thresh)\n');
% load('experiment_05-33-49_07-50-43.mat');
% dists = [];
% for i = 1:length(log_coords1_all)
%     dists = [dists; log_coords1_all{i}(:, 5)'];
% end
% fprintf('\nMeans: %.2f %.2f %2.f %.2f\n', mean(dists));
% fprintf('Median: %.2f %.2f %2.f %.2f\n', median(dists));
% fprintf('Standard deviations: %.2f %.2f %2.f %.2f\n', std(dists));
% fprintf('Max: %.2f %.2f %2.f %.2f\n', max(dists));
% fprintf('Min: %.2f %.2f %2.f %.2f\n\n', min(dists));
% 
% fprintf('\nExperiment 5: (default conv, matlab thresh)\n');
% load('experiment_01-29-57_02-31-58.mat');
% dists = [];
% for i = 1:length(log_coords1_all)
%     dists = [dists; log_coords1_all{i}(:, 5)'];
% end
% fprintf('\nMeans: %.2f %.2f %2.f %.2f\n', mean(dists));
% fprintf('Median: %.2f %.2f %2.f %.2f\n', median(dists));
% fprintf('Standard deviations: %.2f %.2f %2.f %.2f\n', std(dists));
% fprintf('Max: %.2f %.2f %2.f %.2f\n', max(dists));
% fprintf('Min: %.2f %.2f %2.f %.2f\n\n', min(dists));
% 
% 
% 
% fprintf('\nExperiment 6: (default conv, matlab thresh)\n');
% load('experiment_temp_k_fixo.mat');
% dists = [];
% for i = 1:length(log_coords1_all)
%     dists = [dists; log_coords1_all{i}(:, 5)'];
% end
% fprintf('\nMeans: %.2f %.2f %2.f %.2f\n', mean(dists));
% fprintf('Median: %.2f %.2f %2.f %.2f\n', median(dists));
% fprintf('Standard deviations: %.2f %.2f %2.f %.2f\n', std(dists));
% fprintf('Max: %.2f %.2f %2.f %.2f\n', max(dists));
% fprintf('Min: %.2f %.2f %2.f %.2f\n\n', min(dists));
% 
% 
% fprintf('\nExperiment 7: (default conv, matlab thresh)\n');
% load('experiment_temp_k_varia.mat');
% dists = [];
% for i = 1:length(log_coords1_all)
%     dists = [dists; log_coords1_all{i}(:, 5)'];
% end
% fprintf('\nMeans: %.2f %.2f %2.f %.2f\n', mean(dists));
% fprintf('Median: %.2f %.2f %2.f %.2f\n', median(dists));
% fprintf('Standard deviations: %.2f %.2f %2.f %.2f\n', std(dists));
% fprintf('Max: %.2f %.2f %2.f %.2f\n', max(dists));
% fprintf('Min: %.2f %.2f %2.f %.2f\n\n', min(dists));



fprintf('\nExperiment 8: (default conv, matlab thresh)\n');
% load('experiment_K.mat');
load('experiment_K_semantic.mat');
dists = [];
cntExp = 1;
cntImg = 1;
cntTotal = 1;
% distMatrix = zeros(10, 4, 10);
exps = 0.2:0.2:2;
filters = {'Kuwahara', 'Billateral', 'LeastSquares', 'Anisiotropic'};
Image = zeros(length(targets), 1);
K = zeros(length(targets), 1);
Image2 = zeros(length(log_coords1_all), 1);
K2 = zeros(length(log_coords1_all), 1);
Filter = cell(length(targets), 1);
Dist = zeros(length(targets), 1);
Dist_Kuwa = zeros(length(log_coords1_all), 1);
Dist_Bill = zeros(length(log_coords1_all), 1);
Dist_LtSq = zeros(length(log_coords1_all), 1);
Dist_Anis = zeros(length(log_coords1_all), 1);
kTemplate_Direct = zeros(length(targets), 1);
kTemplate_Inverse = zeros(length(targets), 1);
kTarget_Direct = zeros(length(targets), 1);
kTarget_Inverse = zeros(length(targets), 1);
for i = 1:length(log_coords1_all)
    currDists = log_coords1_all{i}(:, 5)';
    kTempDir = log_coords1_all{i}(:, 7)';
    kTempInv = log_coords1_all{i}(:, 8)';
    kTargDir = log_coords1_all{i}(:, 9)';
    kTargInv = log_coords1_all{i}(:, 10)';
    dists = [dists; currDists];
%     distMatrix(cntExp, :, cntImg) = currDists;
    for j = 1:4
        Image(cntTotal) = cntImg;
        K(cntTotal) = exps(cntExp);
        Filter{cntTotal} = filters{j};
        Dist(cntTotal) = currDists(j);
        kTemplate_Direct(cntTotal) = kTempDir(j);
        kTemplate_Inverse(cntTotal) = kTempInv(j);
        kTarget_Direct(cntTotal) = kTargDir(j);
        kTarget_Inverse(cntTotal) = kTargInv(j);
        cntTotal = cntTotal+1;
    end
    Image2(i) = cntImg;
    K2(i) = exps(cntExp);
    Dist_Kuwa(i) = currDists(1);
    Dist_Bill(i) = currDists(2);
    Dist_LtSq(i) = currDists(3);
    Dist_Anis(i) = currDists(4);
    cntExp = cntExp+1;
    if cntExp == 11; cntExp = 1; cntImg = cntImg+1; end
end
T2 = table(Image2, K2, Dist_Kuwa, Dist_Bill, Dist_LtSq, Dist_Anis);
T = table(Image, K, Dist, Filter, kTemplate_Direct, kTemplate_Inverse, kTarget_Direct, kTarget_Inverse);

% exemplos de buscas na tabela
T2(T2.Image2 == 1,:);
T(T.K == 0.2 & strcmp(T.Filter, 'Kuwahara'),:);
find(round(T.Dist, 3) == 12.948);

% plotando figura
im1 = find(round(T.Dist, 4) == 48.6304, 1);
im2 = find(round(T.Dist, 3) == 257.199, 1);
figure; subplot(221); imshow(templates{im1}); title(['Template Image ' num2str(im1)]); subplot(222); imshow(templates{im2}); title(['Template Image ' num2str(im2)]);
subplot(223); imshow(targets{im1}); title(['Target Image ' num2str(im1)]); subplot(224); imshow(targets{im2}); title(['Target Image ' num2str(im2)]);

figure; imshow(targets{im1});

fprintf('\nMeans: %.2f %.2f %2.f %.2f\n', mean(dists));
fprintf('Median: %.2f %.2f %2.f %.2f\n', median(dists));
fprintf('Standard deviations: %.2f %.2f %2.f %.2f\n', std(dists));
fprintf('Max: %.2f %.2f %2.f %.2f\n', max(dists));
fprintf('Min: %.2f %.2f %2.f %.2f\n\n', min(dists));



dbg = 1;