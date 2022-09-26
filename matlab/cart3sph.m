function [r, theta, phi] = cart3sph(q)
x=q(1);
y=q(2);
z=q(3);

[theta,r1] = cart2pol(x,y);
if x==0
    xx=0;
else
    xx=x/cos(theta);
end

[phi,r] = cart2pol(xx,z);
end