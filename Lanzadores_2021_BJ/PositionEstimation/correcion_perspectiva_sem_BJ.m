function [img_out, util_img, utilCropSize]=correcion_perspectiva_sem_BJ(img, yaw, pitch, roll, distz, pre_geo_img, pre_uav_img, perspective_correction_type, printFigs, preprocessing_img)

if perspective_correction_type == 1
    fprintf('\nPerspective Correction with Parametric Logs...\n');
    tic;
    [h, w,~] = size(img);

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
    util_img = img_out;
    utilCropSize = img_out;
    fprintf('Execution time for Perspective Correction: %f s\n', toc);

elseif perspective_correction_type == 2
    
    fprintf('\nPerspective Correction with SURF KeyPoints matching...\n');
    tic;
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
    util_img = img_out;
    utilCropSize = img_out;
    fprintf('Execution time for Perspective Correction: %f s\n', toc);
    
    if printFigs
        figure;
        showMatchedFeatures(pre_geo_img,pre_uav_img,inlierOriginalXY,inlierDistortedXY,'montage')
        figure;
        showMatchedFeatures(pre_uav_img,pre_geo_img,inlierDistortedXY,inlierOriginalXY,'montage') 
        figure;
        imshowpair(pre_geo_img,img_out,'montage')
    end
    
end
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
