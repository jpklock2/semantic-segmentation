function [imout] = useful_area(im ,w, h, theta)
% adapted from:
% https://www.mathworks.com/matlabcentral/fileexchange/48624-rotate-images-with-automatic-cropping
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
xoff = 0;           
yoff = 0;
if ~isnan(delta)
    
    wp =  c * delta * w + s * h;
    hp = -s * delta * w + c * h;
    
elseif ~isnan(gamma)
    
    wp = c * w - s * gamma * h;
    hp = s * w + c * gamma * h;
    
else    
    
    wp = NaN;
    hp = NaN;
    
end

if swapxy
    [wp, hp, xoff, yoff] = deal(hp, wp, yoff, xoff);
end
if tall
    xoff = -xoff;
end

imout = im;
        
[hr, wr, ~] = size(imout);
% dimension are of box on outside of outer pixels
y0 = 1 + ceil((hr - hp)/2 - yoff);
y1 = floor((hr + hp)/2 - yoff);
x0 = 1 + ceil((wr - wp)/2 + xoff);
x1 = floor((wr + wp)/2 + xoff);
imout = imout(y0:y1, x0:x1, :);
end