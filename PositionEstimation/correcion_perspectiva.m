function [distx, disty, Img1]=correcion_perspectiva(img, yaw, pitch, roll, distz)

[h w] = size(img);
% distz=300;
distx =(300/(tand(90-roll)));
disty =(300/(tand(90-pitch)));

% Projection 2D - 3D 
A1=[1 0 -w/2;
    0 1 -h/2;
    0 0    0;
    0 0    1];

% Composed rotation matrix with (RX, RY, RZ)
R_rot = (R_x(roll))*(R_y(pitch))*(R_z(yaw)); %yaw pitch roll

% Translation matrix
Tras=[1 0 0 0;
      0 1 0 0;
      0 0 1 distz;
      0 0 0 1];

% 3D -> 2D 
A2=[650.88   0     w/2   0;
   0     649.43    h/2 0;
   0     0      1     0];

H= A2*(Tras*(R_rot*A1));

tform = projective2d(H');
Img1 = imwarp(img,tform,'cubic');


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

end