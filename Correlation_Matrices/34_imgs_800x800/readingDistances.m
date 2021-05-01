clc;
badAngles = [24 25 32 33 34];

load('templateResults_34_imgs_v3_corrected_useful_area.mat');
allTemplateResults(badAngles, :) = [];
% for i = 1:length(badAngles)
%     allTemplateResults(badAngles(i), 10) = allTemplateResults(badAngles(i), 12);
% end
distsOld = cell2mat(allTemplateResults(:, 12));
distsNew = cell2mat(allTemplateResults(:, 10));
fprintf('\n70%% thresh and 0.3 prob - Old mean dist = %.2f - New mean dist = %.2f\n', mean(distsOld), mean(distsNew));
fprintf('70%% thresh and 0.3 prob - Old std dist =  %.2f - New std dist = %.2f\n', std(distsOld), std(distsNew));

allTemplateResults2 = allTemplateResults;

load('templateResults_34_imgs_v2_04.mat');
allTemplateResults(badAngles, :) = [];
% for i = 1:length(badAngles)
%     allTemplateResults(badAngles(i), 10) = allTemplateResults(badAngles(i), 12);
% end
distsOld = cell2mat(allTemplateResults(:, 12));
distsNew = cell2mat(allTemplateResults(:, 10));
fprintf('\n70%% thresh and 0.4 prob - Old mean dist = %.2f - New mean dist = %.2f\n', mean(distsOld), mean(distsNew));
fprintf('70%% thresh and 0.4 prob - Old std dist =  %.2f - New std dist =  %.2f\n', std(distsOld), std(distsNew));

load('templateResults_34_imgs_v2_75.mat');
% allTemplateResults(badAngles, :) = [];
distsOld = cell2mat(allTemplateResults(:, 12));
distsNew = cell2mat(allTemplateResults(:, 10));
fprintf('\n75%% thresh and 0.3 prob - Old mean dist = %.2f - New mean dist = %.2f\n', mean(distsOld), mean(distsNew));
% fprintf('75%% thresh and 0.3 prob - Old std dist =  %.2f - New std dist =  %.2f\n', std(distsOld), std(distsNew));

load('templateResults_34_imgs_v2_04_75.mat');
% allTemplateResults(badAngles, :) = [];
distsOld = cell2mat(allTemplateResults(:, 12));
distsNew = cell2mat(allTemplateResults(:, 10));
fprintf('\n75%% thresh and 0.4 prob - Old mean dist = %.2f - New mean dist = %.2f\n\n', mean(distsOld), mean(distsNew));

load('templateResults_34_imgs_v2.mat');
distsOld = cell2mat(allTemplateResults(badAngles, 12));
distsNew = cell2mat(allTemplateResults(badAngles, 10));
fprintf('\nOld useful area, images 24, 25, 32, 33, 34 - Old mean dist = %.2f - New mean dist = %.2f\n', mean(distsOld), mean(distsNew));

load('templateResults_34_imgs_v3_corrected_useful_area.mat');
distsOld = cell2mat(allTemplateResults(badAngles, 12));
distsNew = cell2mat(allTemplateResults(badAngles, 10));
fprintf('\nNew useful area, images 24, 25, 32, 33, 34 - Old mean dist = %.2f - New mean dist = %.2f\n\n', mean(distsOld), mean(distsNew));

% fprintf('\n70%% thresh and 0.3 prob - Old mean dist = %.2f - New mean dist = %.2f\n', mean(log_coords1(:, 8)), mean(log_coords1(:, 5)));
% fprintf('70%% thresh and 0.3 prob - Old std dist =  %.2f - New std dist = %.2f\n\n', std(log_coords1(:, 8)), std(log_coords1(:, 5)));

