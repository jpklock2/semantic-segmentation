function [mask] = removeResiduals(mask, classes)

    % removing residuals based on shortest distance
    mask = round(mask);
    change = mask(~ismember(mask, [0; unique(classes)]));
    mantain = unique(classes);
    mantainMatrix = ones(size(change, 1), size(mantain, 1)).*mantain';
    dists = sqrt((change - mantainMatrix).^2);
    [~, min_dists] = max((dists == min(dists, [], 2)), [], 2);
    replace = mantainMatrix((1:size(dists, 2)) == min_dists);
    mask(~ismember(mask, [0; unique(classes)])) = replace;

end

