function [img_out, util_img, Lmask, util_Lmask, mask, utilCropSize, currImg, currPlotImage, currSegmentation]=correcion_perspectiva_segmentation(img, yaw, pitch, roll, distz, Lmask, mask, currImg, currPlotImage, currSegmentation,pre_geo_img, pre_uav_img, perspective_correction_type,preprocessing_img )
if perspective_correction_type == 1
    
[h, w,~] = size(img);
oldImg = zeros(size(img));

distx =(300/(tand(90-roll)));
disty =(300/(tand(90-pitch)));

% Projection 2D - 3D 
A1=[1 0 -w/2;
    0 1 -h/2;
    0 0    1;
    0 0    1];

% Composed rotation matrix with (RX, RY, RZ)
R_rot = (R_x(roll))*(R_y(pitch))*(R_z(yaw)); %yaw pitch roll

% Translation matrix
Tras=[1 0 0 0;
      0 1 0 0;
      0 0 1 distz;
      0 0 0 1];

% 3D -> 2D 
A2=[distz   0     w/2   0;
   0     distz    h/2 0;
   0     0      1     0];

H = A2*(Tras*(R_rot*A1));
tform = projective2d(H');
if preprocessing_img
    img_out = imwarp(pre_uav_img,tform,'bilinear');
    img_out_copy = imwarp(img,tform,'bilinear');
else
    img_out = imwarp(img,tform,'bilinear');
    img_out_copy = img_out;
end
if exist('Lmask','var')
    Lmask = imwarp(Lmask,tform,'nearest');
    mask = imwarp(mask,tform, 'nearest');
    currImg = imwarp(currImg,tform,'bilinear');
    currSegmentation = imwarp(currSegmentation,tform,'bilinear');
    currPlotImage = imwarp(currPlotImage,tform,'bilinear');
end


    
elseif perspective_correction_type == 2
    ptsOriginalSURF  = detectSURFFeatures(pre_geo_img,'MetricThreshold',400,'NumOctaves', 16 ,'NumScaleLevels',20);
    ptsDistortedSURF = detectSURFFeatures(pre_uav_img,'MetricThreshold',400,'NumOctaves', 16,'NumScaleLevels',20);

    [featuresOriginalSURF,  validPtsOriginalSURF]  = extractFeatures(pre_geo_img,  ptsOriginalSURF,'Method','SURF');
    [featuresDistortedSURF, validPtsDistortedSURF] = extractFeatures(pre_uav_img, ptsDistortedSURF,'Method','SURF');
 
    indexPairsSURF = matchFeatures(featuresOriginalSURF, featuresDistortedSURF,'MatchThreshold',10, 'MaxRatio', 0.99,'Metric','SAD');

    matchedOriginalSURF  = validPtsOriginalSURF(indexPairsSURF(:,1),:);
    matchedDistortedSURF = validPtsDistortedSURF(indexPairsSURF(:,2),:);
  
    [tform,inlierDistortedXY,inlierOriginalXY] = estimateGeometricTransform(matchedDistortedSURF, ...
        matchedOriginalSURF ,'similarity','MaxDistance',15,'Confidence',90,'MaxNumTrials',5000);%0 17 en maxdistance con 1500 es una buena metrica(17000 por defecto)
     
    puntos_VANT=inlierDistortedXY;
    puntos_georef=inlierOriginalXY;
    
    outputView = imref2d(size(pre_geo_img));
    if preprocessing_img 
        img_out = imwarp(pre_uav_img,tform,'bilinear');
        img_out_copy = imwarp(img,tform,'bilinear');
    else
        img_out  = imwarp(img,tform,'bilinear');
        img_out_copy = img_out;
    end
end
    
theta = atan2(tform.T(2,1), tform.T(2,2)) * 180 / pi;

[y0, x0, y1, x1] = useful_area_v2(oldImg, currImg, theta, 1, 0);
[y00, x00, y11, x11] = useful_area_v2(oldImg, currImg, theta, 0, 0);
area0 = (y1-y0)*(x1-x0);
area1 = (y11-y00)*(x11-x00);
if (abs(theta) > 45 && abs(theta) < 135) || (abs(theta) > 225 && abs(theta) < 315)
    y0 = y00; y1 = y11; x0 = x00; x1 = x11;
end

util_img = img_out_copy;

[wp, hp, xoff, yoff] = useful_area(w, h, yaw);
[hr, wr, ~] = size(util_img);
y0old = 1 + ceil((hr - hp)/2 - yoff);
y1old = floor((hr + hp)/2 - yoff);
x0old = 1 + ceil((wr - wp)/2 + xoff);
x1old = floor((wr + wp)/2 + xoff);
xcrdsOld = [x0old x0old x1old x1old x0old];
ycrdsOld = [y0old y1old y1old y0old y0old];
% plot(xcrdsOld, ycrdsOld, 'g-');

[wp, hp, xoff, yoff] = useful_area(w, h, theta);
[hr, wr, ~] = size(util_img);
y0old = 1 + ceil((hr - hp)/2 - yoff);
y1old = floor((hr + hp)/2 - yoff);
x0old = 1 + ceil((wr - wp)/2 + xoff);
x1old = floor((wr + wp)/2 + xoff);
xcrdsOld = [x0old x0old x1old x1old x0old];
ycrdsOld = [y0old y1old y1old y0old y0old];
% plot(xcrdsOld, ycrdsOld, 'b-');

utilCropSize = [y0 x0 y1 x1];
util_img = img_out(y0:y1, x0:x1, :);
if exist('Lmask','var')
    util_Lmask = Lmask(y0:y1, x0:x1, :);
%     util_mask = mask(y0:y1, x0:x1, :);
end

%% Matrix for Yaw-rotation about the Z-axis
function [R] = R_z(psi)  
    R = [cosd(psi)   -sind(psi) 0  0;
         sind(psi)  cosd(psi)   0  0;
         0           0          1  0;
         0           0          0  1];
end

%% Matrix for Pitch-rotation about the Y-axis
function [R] = R_y(theta) 
    R = [cosd(theta)    0  -sind(theta)  0;
         0              1   0            0;
         sind(theta)    0   cosd(theta)  0;
         0              0       0        1];
end
%% Matrix for Roll-rotation about the X-axis
function [R] = R_x(phi) 
    R = [1  0           0          0;
         0  cosd(phi)   -sind(phi) 0;
         0  sind(phi)  cosd(phi)   0;
         0     0           0       1];
end

%%
function [wp, hp, xoff, yoff] = useful_area(w, h, theta)
% adapted from https://www.mathworks.com/matlabcentral/fileexchange/48624-rotate-images-with-automatic-cropping
a = w / h;              
tall = w < h;
if tall
    [w, h] = deal(h, w);    % invert aspect ratio if necessary
end
csgn = cosd(theta);         % save original rotation
ssgn = sind(theta);         % angle
c = abs(csgn);              % symmetries allow most of computation to use
s = abs(ssgn);              % absolute rotation angles
bigrot = s > c;
if bigrot                   % take complement of angle if necessary
    [s, c, ssgn, csgn] = deal(c, s, csgn, ssgn);
end
swapxy = xor(tall, bigrot);
% Now working with c, s in range 0 to 1/sqrt(2) and r in range 0 to 1
r = h / w;                  % inverse original aspect ratio
sin2th = 2 * c * s;
cos2th = c*c - s*s;
gamma = NaN;
delta = NaN;
fourpt = false;
 
if sin2th < r           % four contact is max area solution?
    fourpt = true;
    if cos2th > eps     % test for special 45 degree, r=1 case
        delta = (1 - r*sin2th)/cos2th; % general 4-contact case
    else
        delta = 0;      % special 4-contact case
    end
else
    delta = r * cos2th / sin2th; % 2-contact max area solution
end

% use delta or gamma to find output width and height
xoff = 0;           % default position offset
yoff = 0;
if ~isnan(delta)
    
    wp =  c * delta * w + s * h;
    hp = -s * delta * w + c * h;
    
elseif ~isnan(gamma)
    
    wp = c * w - s * gamma * h;
    hp = s * w + c * gamma * h;
    
else    % four point solution requested, but not possible
    
    wp = NaN;
    hp = NaN;
    
end
% reverse transform carried out at start
if swapxy
    [wp, hp, xoff, yoff] = deal(hp, wp, yoff, xoff);
end
if tall
    xoff = -xoff;
end

end

%%
function [y0, x0, y1, x1] = useful_area_v2(img, currImg, theta, mantainAspectRatio, plotFig)

[h, w,~] = size(img);
    
if plotFig
figure; imshow(currImg);
hold on;
end
boxin.width = size(currImg, 2);
boxin.height = size(currImg, 1);
boxout = imRotateCrop(boxin, theta, 'AspectRatio', 'maxArea', 'Position', 0);
xcentre = (1+size(currImg, 2))/2 + boxout.xshift;
ycentre = (1+size(currImg, 1))/2 + boxout.yshift;
w2 = boxout.width/2;
h2 = boxout.height/2;
xcrds = xcentre + [w2 w2 -w2 -w2 w2];   % box corner coords
ycrds = ycentre + [h2 -h2 -h2 h2 h2];
upperLim = floor(ycentre+h2);
bottomLim = floor(ycentre-h2);
if plotFig
hold on;
plot(xcrds, ycrds, 'g-', 'LineWidth', 3);
end

[h1, w1, ~] = size(currImg);
w2 = w1/2;
h2 = h1/2;

dg = sqrt(w2^2 + h2^2);
angInc = asin(h2/dg)*180/pi;
xc = (1+w1)/2;
yc = (1+h1)/2;

if mantainAspectRatio
    diffH = h1 - h;
    diffW = w1 - w;
else
    diffH = 0;
    diffW = 0;
end

c1 = polyfit([xc w1-diffW], [yc 1+diffH], 1);
yc1 = round(linspace(xc, w1, floor(0.95*dg)));
yc1 = c1(1)*yc1 + c1(2);

c2 = polyfit([xc 1+diffW], [yc 1+diffH], 1);
yc2 = round(linspace(xc, 1, floor(0.95*dg)));
yc2 = c2(1)*yc2 + c2(2);

hi = max([round(yc1(end)) 1]);

c3 = polyfit([xc 1+diffW], [yc h1-diffH], 1);
yc3 = round(linspace(xc, 1, floor(0.99*dg)));
yc3 = c3(1)*yc3 + c3(2);

c4 = polyfit([xc w1-diffW], [yc h1-diffH], 1);
yc4 = round(linspace(xc, w1, floor(0.99*dg)));
yc4 = c4(1)*yc4 + c4(2);

hf = min([round(yc3(end)) h1]);

diffW = 0;

if plotFig
l1 = line([xc w1-diffW], [yc hi],'Color','red','LineWidth', 3);
l2 = line([xc 1+diffW], [yc hi],'Color','red','LineWidth', 3);
l3 = line([xc 1+diffW], [yc hf],'Color','red','LineWidth', 3);
l4 = line([xc w1-diffW], [yc hf],'Color','red','LineWidth', 3);
end

xd1 = round(linspace(xc, w1-diffW, floor(0.99*dg)));
yd1 = round(linspace(yc, hi, floor(0.99*dg)));

xd2 = round(linspace(xc, 1+diffW, floor(0.99*dg)));
yd2 = round(linspace(yc, hi, floor(0.99*dg)));

xd3 = round(linspace(xc, 1+diffW, floor(0.99*dg)));
yd3 = round(linspace(yc, hf, floor(0.99*dg)));

xd4 = round(linspace(xc, w1-diffW, floor(0.99*dg)));
yd4 = round(linspace(yc, hf, floor(0.99*dg)));

flagQ1 = 0; flagQ2 = 0; flagQ3 = 0; flagQ4 = 0;
for i = 1:length(xd1)
    if all(currImg(yd1(i), xd1(i), :) == 0) && ~flagQ1 %&& (yd1(i) < yc+h2)
        q1 = [xd1(i) yd1(i)];
        flagQ1 = 1;
    end
    
    if all(currImg(yd2(i), xd2(i), :) == 0) && ~flagQ2 %&& (yd2(i) < yc+h2)
        q2 = [xd2(i) yd2(i)];
        flagQ2 = 1;
    end
    
    if all(currImg(yd3(i), xd3(i), :) == 0) && ~flagQ3 %&& (yd3(i) < yc+h2)
        q3 = [xd3(i) yd3(i)];
        flagQ3 = 1;
    end
    
    if all(currImg(yd4(i), xd4(i), :) == 0) && ~flagQ4 %&& (yd4(i) < yc+h2)
        q4 = [xd4(i) yd4(i)];
        flagQ4 = 1;
    end
end

if flagQ1 == 0
    q1 = [w1-diffW hi];
end

if flagQ2 == 0
    q2 = [1+diffW hi];
end

if flagQ3 == 0
    q3 = [1+diffW hf];
end

if flagQ4 == 0
    q4 = [w1-diffW hf];
end

y0 = max([q1(2) q2(2) bottomLim]);
y1 = min([q3(2) q4(2) upperLim]);
x0 = max(q2(1), q3(1));
x1 = min(q1(1), q4(1));

if plotFig
xcrds = [x0 x0 x1 x1 x0];
ycrds = [y0 y1 y1 y0 y0];
plot(xcrds, ycrds, 'b-', 'LineWidth', 3);
end

end

end