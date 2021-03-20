if (useCropped)
    load(['./Images/Test/Cropped/Classes/cropped_im',num2str(m),'.mat']);
else
    load(['./Images/Test/Original_Dev/Classes/im',num2str(m),'.mat']);
end

if increaseSp
    Lmask = zeros(size(L2));
    for i = 1:length(myClasses)
        currC = myClasses{i};
        Lmask(ismember(L2, currC)) = i;
    end
    myClassesNew = cell(size(myClasses));
    for i = 1:N
        [GC,GR] = groupcounts(Lmask(L == i));
        class = GR(GC == max(GC));
        myClassesNew{class(1)} = [myClassesNew{class(1)}; i];
    end
    myClasses = myClassesNew;
else
    L = L2;
    N = N2;
end

if isempty(myClasses)
    mapImage = rgbImage;
    classesImage = tempImage;
    myClasses = {};
    for class = 1:length(myColors)
        spClasses = [];
        for sp = 1:N
            [xL, yL] = find(L == sp);
            currClass = classesImage(xL(1), yL(1), :);
            currClass = reshape(currClass, [1 3]);
            if all(currClass == myColors(class, :))
                spNumber = sp;
                spClasses = [spClasses; spNumber];
            end
        end
        myClasses{class} = spClasses;
    end
%     save myClasses.mat mapImage classesImage L2 N2 myClasses
end

if combine
    myClasses{1} = [myClasses{1}; myClasses{2}; myClasses{3}];
    myClasses{2} = [myClasses{4}; myClasses{5}; myClasses{9}; myClasses{12}];
    myClasses{3} = [myClasses{7}; myClasses{8}; myClasses{11}];
    myClasses{4} = [myClasses{6}; myClasses{10}; myClasses{13}; myClasses{14}];
    myClasses = myClasses(1:4);
end

idx = label2idx(L);
BW = boundarymask(L);
allC = 1:N;
classes = zeros(N, 1);
classesTemp = zeros(N, length(myClasses));
for i = 1:length(myClasses)
    currC = myClasses{i};
    classesTemp(currC, i) = 1;
    classes(currC) = i;
    allC(currC) = 0;
end
disp(allC(allC ~= 0));
dbg = 1;
% U = classesTemp;
% classes = classes';