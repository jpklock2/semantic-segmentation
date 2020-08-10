function [img, util_img]=correcion_perspectiva(img, yaw, pitch, roll, distz)

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
img = imwarp(img,tform,'bilinear');

[wp, hp, xoff, yoff] = useful_area(w, h, yaw);
util_img = img;
[hr, wr, ~] = size(util_img);

y0 = 1 + ceil((hr - hp)/2 - yoff);
y1 = floor((hr + hp)/2 - yoff);
x0 = 1 + ceil((wr - wp)/2 + xoff);
x1 = floor((wr + wp)/2 + xoff);

util_img = util_img(y0:y1, x0:x1, :);

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

end