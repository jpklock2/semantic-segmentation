function [fig] = plotMask(mask, classes, rec, adjacencies, indexes, maskC, crop, drone)

if ~exist('maskC','var')
    maskC = zeros(size(mask));
    for mx = 1:length(classes)
        maskC(mask == mx) = classes(mx);
    end
end

BW = boundarymask(mask);
img = imoverlay(ind2rgb(im2uint8(mat2gray(maskC)), parula(256)), BW, 'k');

if exist('adjacencies','var') && exist('indexes','var') && ~isempty(adjacencies)
lines = [{}];
for i = 1:length(adjacencies)
    img = insertText(img, [adjacencies{i}.centX adjacencies{i}.centY], num2str(i), 'FontSize', 18, 'TextColor', 'r', 'BoxOpacity', 0);
    if ~isempty(adjacencies{i}.adj)
        for j = 1:length(adjacencies{i}.adj)
            currAdj = adjacencies{i}.adj(j);
            if classes(indexes(i)) ~= classes(indexes(currAdj))
                lines = [lines; [adjacencies{i}.centY adjacencies{i}.centX adjacencies{currAdj}.centY adjacencies{currAdj}.centX]];
            end
        end
    end
end

for i = 1:length(lines)
    line([lines{i}(2) lines{i}(4)], [lines{i}(1) lines{i}(3)], 'Color', 'r', 'Marker', 'o');
end
end

fig = figure; imshow(img);
hold on;
rectangle('Position',[rec(2), rec(1), rec(4)-rec(2), rec(3)-rec(1)],'LineWidth',2,'EdgeColor','r');

if exist('crop','var')
    rectangle('Position',[crop(2), crop(1), crop(4)-crop(2), crop(3)-crop(1)],'LineWidth',2,'EdgeColor','k');
end

if exist('drone','var')
    rectangle('Position',[drone(2), drone(1), 2, 2],'LineWidth',5,'EdgeColor','w','Curvature',[1 1]);
end
    
dbg = 1;

end
