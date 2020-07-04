% cria minha MF customizada

numInp = size(pixelsOri, 2);
numOutp = size(classesOri, 2);
cluster_n = KF;
U_temp = U';

fismat = sugfis('Name','sugenoClasses');
fismat.DisableStructuralChecks = true;
minX = zeros(1, numInp); % min(pixelsOri);
maxX = ones(1, numInp); % max(pixelsOri);
ranges = fuzzy.internal.utility.updateRangeIfEqual([minX ; maxX]');
mfTemplate = fismf;
mfTemplate.Type = 'gaussmf';
mf = repmat(mfTemplate,[1 cluster_n]);
var = repmat(fisvar,[1 numInp]);
for i = 1:numInp
    inputName = ['in' num2str(i)];
    var(i).Name = inputName;
    var(i).Range = ranges(i,:);

    % Loop through and add mf's
    for j = 1:1:cluster_n
        mf(j).Name = ['in' num2str(i) 'cluster' num2str(j)];
        mf(j).Parameters = computemfparams ('gaussmf', pixelsOri(:,i), U_temp(j,:)', centroidsFinal(j,i));
    end   
    var(i).MembershipFunctions = mf;
end
fismat.Inputs = var;

minX = min(classesOri);
maxX = max(classesOri);
ranges = fuzzy.internal.utility.updateRangeIfEqual([minX ; maxX]');
var = repmat(fisvar,[1 numOutp]);

mfS = repmat(fismf('linear',zeros(1,length(fismat.Inputs)+1)),[1 cluster_n]);
% Loop through and add outputs
for i=1:1:numOutp
    outputName = ['out' num2str(i)];
    var(i).Name = outputName;
    var(i).Range = ranges(i,:);

    % Loop through and add mf's
    for j = 1:1:cluster_n
        mfS(j).Name = ['out' num2str(i) 'cluster' num2str(j)];
        mfS(j).Parameters = computemfparams ('linear', [pixelsOri classesOri(:,i)]);
    end
    var(i).MembershipFunctions = mfS;
end
fismat.Outputs = var;

% Create rules
ruleList = ones(cluster_n, numInp+numOutp+2);
for i = 2:1:cluster_n
    ruleList(i,1:numInp+numOutp) = i;    
end
fismat = addRule(fismat,ruleList);

fismat = enableCheckWithoutConstruction(fismat);
inFis = fismat;

function mfparams = computemfparams(mf,x,m,c)

switch lower(mf)
    
    case 'gaussmf'
        sigma = invgaussmf4sigma (x, m, c);
        mfparams = [sigma, c];
        
    case 'linear'
        [N, dims] = size(x);
        xin = [x(:,1:dims-1) ones(N,1)];
        xout = x(:, dims);
        b = xin \ xout;
        mfparams = b';
        
    % Don't need 'otherwise' case since 'mf' is hard coded to 'gaussmf'.
end

end
