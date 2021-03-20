function [mask] = removeResidualsMask(mask)

    % removing residuals based on shortest distance
%     mask = round(mask);
    [y, x] = find(mod(mask, 1) ~= 0);
    for i = 1:length(x)
        currMask = mask(max(1, y(i)-2):min(size(mask, 1), y(i)+2), max(1, x(i)-2):min(size(mask, 2), x(i)+2));
        change = [0; unique(currMask(mod(currMask, 1) == 0))];
        dists = sqrt((mask(y(i), x(i)) - change).^2);
        minIdx = find(dists == min(dists), 1);
        mask(y(i), x(i)) = change(minIdx);
    end
end

