function [classesTempJoin, KcJoin, UcJoin, centroidsJoin] = joinClasses(tree, Kc, classesTemp, centroids, Uc, pixels, thresh)

Z = tree;
Z = transz(Z);
thresh = 0.33*max([Z(:, 3)]);
% thresh = 0.33/max([Z(:, 3); 1]);
usable = Z(Z(:, 3) <= thresh, 1:2);
% currGroups = 1:Kc;

groups = [{}];
for i = size(usable, 1):-1:1
    if isempty(groups)
        groups = [groups; usable(i, :)];
%         currGroups(ismember(currGroups, usable(i, :))) = [];
    else
        newGroup = 1;
        for j = 1:length(groups)
           if  sum(ismember(groups{j}, usable(i, :))) > 0
               groups{j} = unique([groups{j}(:); usable(i, :)']);
%                currGroups(ismember(currGroups, groups{j})) = [];
               newGroup = 0;
           end
        end
        if newGroup == 1
            groups = [groups; usable(i, :)];
        end
    end
end

% KcJoin = length(currGroups);
classesTempJoin = classesTemp;
allIdx = 1:Kc;
eliminate = [];
for i = 1:length(groups)
    classesTempJoin(ismember(classesTempJoin, groups{i})) = groups{i}(1);
    Uc(:, groups{i}(1)) = sum(Uc(:, groups{i}), 2);
    eliminate = [eliminate; groups{i}(2:end)];
end
if ~isempty(eliminate)
    allIdx(eliminate) = [];
    UcJoin = Uc(:, allIdx);
    centroidsJoin = ((pixels'*UcJoin)./sum(UcJoin))';
else
    UcJoin = Uc;
    centroidsJoin = centroids;
end

fixList = unique(classesTempJoin);
KcJoin = length(fixList);
classesTempJoin = classesTempJoin + 2*Kc;
for i = 1:length(fixList)
    classesTempJoin(ismember(classesTempJoin, fixList(i)+2*Kc)) = i;
end

dbg = 1;

% numLeaves = length(leafOrder);
% [X, Y] = orderTree(numLeaves, Z, leafOrder, 1);
% [A, B] = getCoordinatesAndColors(numLeaves, X, Y, Z, theGroups, groups, cmap);

end

% Tree functions
function Z = transz(Z)
    %TRANSZ Translate output of LINKAGE into another format.
    %   This is a helper function used by DENDROGRAM and COPHENET.
    %   For each node currently labeled numLeaves+k, replace its index by
    %   min(i,j) where i and j are the nodes under node numLeaves+k.

    %   In LINKAGE, when a new cluster is formed from cluster i & j, it is
    %   easier for the latter computation to name the newly formed cluster
    %   min(i,j). However, this definition makes it hard to understand
    %   the linkage information. We choose to give the newly formed
    %   cluster a cluster index M+k, where M is the number of original
    %   observation, and k means that this new cluster is the kth cluster
    %   to be formed. This helper function converts the M+k indexing into
    %   min(i,j) indexing.

    numLeaves = size(Z, 1) + 1;

    for i = 1:(numLeaves - 1)
        if Z(i, 1) > numLeaves
            Z(i, 1) = traceback(Z, Z(i, 1));
        end

        if Z(i, 2) > numLeaves
            Z(i, 2) = traceback(Z, Z(i, 2));
        end

        if Z(i, 1) > Z(i, 2)
            Z(i, 1:2) = Z(i, [2 1]);
        end
    end
end

function a = traceback(Z, b)
    numLeaves = size(Z, 1) + 1;

    if Z(b - numLeaves, 1) > numLeaves
        a = traceback(Z, Z(b - numLeaves, 1));
    else
        a = Z(b - numLeaves, 1);
    end

    if Z(b - numLeaves, 2) > numLeaves
        c = traceback(Z, Z(b - numLeaves, 2));
    else
        c = Z(b - numLeaves, 2);
    end

    a = min(a, c);
end

function [X, Y, perm] = orderTree(numLeaves, Z, leafOrder, check)
    % Initializes X such that there will be no crossing
    Y = zeros(numLeaves, 1);

    if isempty(leafOrder)
        r = Y;
        W = arrangeZIntoW(numLeaves, Z);
        [X, perm] = fillXFromW(numLeaves, W, r);
    else % if a leaf order is specified
        X(leafOrder) = 1:numLeaves; % get X based on the specified order
        if (check) % check whether leafOrder will have crossing branch
            checkCrossing(Z(:, 1:2), leafOrder);
        end
        perm = leafOrder;
    end
end

function [X, perm] = fillXFromW(numLeaves, W, r)
    X = 1:numLeaves; %the initial points for observation 1:n
    g = 1;

    for k = 1:numLeaves - 1
        i = W(k, 1); % the left node in W(k,:)

        if ~r(i)
            X(i) = g;
            g = g + 1;
            r(i) = 1;
        end

        i = W(k, 2); % the right node in W(k,:)

        if ~r(i)
            X(i) = g;
            g = g + 1;
            r(i) = 1;
        end
    end
    perm(X) = 1:numLeaves;
end

function [A, B, col, X, Y] = getCoordinatesAndColors(numLeaves, X, Y, Z)
    A = zeros(4, numLeaves - 1);
    B = A;
    col = zeros(numLeaves - 1, 3);

    for n = 1:(numLeaves - 1)
        i = Z(n, 1); j = Z(n, 2); w = Z(n, 3);
        A(:, n) = [X(i) X(i) X(j) X(j)]';
        B(:, n) = [Y(i) w w Y(j)]';
        X(i) = (X(i) + X(j)) / 2;
        Y(i) = w;

    end
end

function checkCrossing(Z, order)
    % check whether the give Tree will have crossing branches
    % with the given permutation vector

    numBranches = size(Z, 1);
    numLeaves = numBranches + 1;
    % reorder the tree
    perm = order(:);
    % XPos is the position indices for leaves 1:numLeaves.
    XPos(perm) = 1:numLeaves;
    Z0 = Z; % keep the original tree
    % renumber the leave nodes in Z such that number N represents the
    % Nth nodes in the plot
    Z = XPos(Z);
    % Check if the reordered tree structure leads to a
    % tree with no crossing branches
    minPos = 1:numLeaves;
    maxPos = 1:numLeaves;
    sz = ones(numLeaves, 1);

    for i = 1:numBranches
        currentMinPos = min(minPos(Z(i, :)));
        currentMaxPos = max(maxPos(Z(i, :)));
        currentSize = sum(sz(Z(i, :)));

        if currentMaxPos - currentMinPos + 1 ~= currentSize
            warning(message('stats:dendrogram:CrossingBranches'));
            break;
        end

        j = XPos(Z0(i, 1)); % j is the cluster number for the newly formed cluster.
        % Note that we can't use j = XPos(min(Z0(i,:))),
        % Because when not all of the points are shown, the value in the
        % first column of Z may be bigger than the value in the second column.
        minPos(j) = currentMinPos;
        maxPos(j) = currentMaxPos;
        sz(j) = currentSize;
    end
end
