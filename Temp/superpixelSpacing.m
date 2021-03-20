function [L, k, nodeRows, nodeCols] = superpixelSpacing(im, k)

% Itest = imread('../original/IMG_4563.JPG');
% Itest = imresize(Itest,0.25,'box');

% s1 = size(Itest, 1);
% s2 = size(Itest, 2);
% K = 1000;
% gridX = round(s1/sqrt(s1*s2/K));
% gridY = round(s2/sqrt(s1*s2/K));
% K = gridX * gridY;

tic;

% k = 1000;
% im = rgb2lab(Itest); 
[rows, cols, chan] = size(im);
% Nominal spacing between grid elements assuming hexagonal grid
S = sqrt(rows*cols / (k * sqrt(3)/2));
% Get nodes per row allowing a half column margin at one end that alternates
% from row to row
nodeCols = round(cols/S - 0.5);
% Given an integer number of nodes per row recompute S
S = cols/(nodeCols + 0.5); 
% Get number of rows of nodes allowing 0.5 row margin top and bottom
nodeRows = round(rows/(sqrt(3)/2*S));
vSpacing = rows/nodeRows;
% Recompute k
k = nodeRows * nodeCols;
% Allocate memory and initialise clusters, labels and distances.
C = zeros(2,k);          % Cluster centre data  1:3 is mean Lab value,
                         % 4:5 is row, col of centre, 6 is No of pixels
% Initialise clusters on a hexagonal grid
kk = 1;
r = vSpacing/2;
for ri = 1:nodeRows
    % Following code alternates the starting column for each row of grid
    % points to obtain a hexagonal pattern. Note S and vSpacing are kept
    % as doubles to prevent errors accumulating across the grid.
    if mod(ri,2), c = S/2; else, c = S;  end
    for ci = 1:nodeCols
        cc = round(c); rr = round(r);
        C(1:2, kk) = [cc; rr];
        c = c+S;
        kk = kk+1;
    end
    r = r+vSpacing;
end
% Now perform the clustering.  10 iterations is suggested but I suspect n
% could be as small as 2 or even 1
S = round(S);  % We need S to be an integer from now on

flag = 1;
cnt = 0;
while flag
    
[L,N] = superpixels(im, 2*k+round(0.5*k*cnt));
% figure; imshow(imoverlay(Itest, boundarymask(L),'cyan'));
Cy = C(1,:);
Cx = C(2,:);
sizes = [Cx; Cy];
cls = zeros(1, N);
for i = 1:N
    [x, y] = find(L == i);
    x_med = round((max(x)+min(x))/2);
    y_med = round((max(y)+min(y))/2);
    dists = dist([x_med y_med], sizes);
    currCls = find(dists == min(dists));
    cls(i) = currCls(1);
end
cnt = cnt+1;

if (size(unique(cls), 2) == k)
    flag = 0;
end

end

idx = label2idx(L);
for i = 1:N
    currIdx = idx{i};
    L(currIdx) = cls(i);
end

toc
figure; imshow(imoverlay(im, boundarymask(L),'cyan'));
for i = 1:size(C, 2)
    text(C(1,i), C(2,i), num2str(i), 'FontSize', 6, 'Color', 'red', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    pause(0.05);
end

% Nvec = [];
% comp = 10;
% for j = 5:50
%     [L,N] = superpixels(Itest, 2000, 'Method', 'slic', 'Compactness', j);
%     Nvec = [Nvec; N];
%     if (N == K)
%         comp = j;
%         break;
%     end
% %     figure; imshow(imoverlay(Itest, boundarymask(L),'cyan'));
% end
% 
% for i = 1:N
%     [x, y] = find(L == i);
%     x_med = round((max(x)+min(x))/2);
%     y_med = round((max(y)+min(y))/2);
%     text(y_med, x_med, num2str(i), 'FontSize', 6, 'Color', 'red', 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
%     pause(0.01);
% end

end