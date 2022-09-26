function F = Flambda(q,lambda)
[theta, r] = cart2pol(q(1),q(2));
F = [r*cos(lambda*theta), r*sin(lambda*theta)];
F = F.';
end