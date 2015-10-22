function [angle,displ,midpoint,seglen] = f_gt(x1,y1,x2,y2)

midpoint = round([x1+x2 y1+y2]/2);

p = [x2-x1 y2-y1];
seglen = norm(p);
p = p/seglen;
angle = cart2pol(p(1),p(2));
% angle = ang(p(1),p(2)); % in [0, 2*pi)
if angle > pi
    angle = angle-pi;
end
p = 0.5*[x1+x2 y1+y2];
displ = dot(p,[-sin(angle) cos(angle)]);
if angle > pi/2
    displ = -displ;
end

end