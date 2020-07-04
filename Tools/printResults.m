clc;
% alfas = [1e-2; 1e-3; 1e-4; 1e-5; 1e-6; 1e-7];
% eps = [100; 200; 500];
% alfas = [1e-1; 0.05; 1e-2; 1e-3; 1e-4; 1e-5; 1e-6];
alfas = 1:10;
eps = [100];
% load accTest7.mat

% load acc_F0_CS0_FE1_C1_A0_CM0.mat
% load acc_F0_CS0_FE1_C3_A0_CM0.mat
% load acc_F0_CS0_FE1_C1_A0_CM1.mat
% load acc_F0_CS0_FE1_C3_A0_CM1.mat

% load acc_F0_CS0_FE1_C1_A1_CM0.mat
% load acc_F0_CS0_FE1_C3_A1_CM0.mat
% load acc_F0_CS0_FE1_C1_A1_CM1.mat
% load acc_F0_CS0_FE1_C3_A1_CM1.mat

% load crossV.mat
% load accCropped.mat
% load accOriginal.mat
% load accOriginalPreprocess1.mat
load testGLCM_old.mat


writematrix(100*accCombinedAll, 'pca_2.xlsx')


% myAccss = myAccs;
% myAccss(:, 1:2) = myAccss(:, 1:2).*100;
myAccss = 100.*accCombinedAll;
% myAccss = 100.*myAccss;
myAccss = round(myAccss, 2);

for i=1:size(myAccss, 2)
    fprintf(['\nCOLUNA ' num2str(i) '\n']);
    for j = 1:size(myAccss, 1)
        myStr = [num2str(myAccss(j,i)) '\n'];
        myStr(strfind(myStr, '.')) = ',';
        fprintf(myStr);
    end
end