function dF = dFlambda(q,h,lambda)
[theta, r] = cart2pol(q(1),q(2));


dCartPol = zeros(2,2);
dCartPol(1,:) = 1/r * [r*cos(theta), r*sin(theta)];
dCartPol(2,:) = 1/r * [-sin(theta), cos(theta)];

h2 = dCartPol*h;
h2(2) = lambda * h2(2);

dPolCart = zeros(2,2);
dPolCart(1,:) = [cos(theta), -r*sin(theta)];
dPolCart(2,:) = [sin(theta), r*cos(theta)];

dF = dPolCart*h2;
end