function [ct T] = shooting_geod(c,h,T,t,a)

q = plf_to_srvf(c,T);
k = deriv_srvf(c,h,T);
[d n] = size(q);

Fq = zeros(d,n);
dFk = zeros(d,n);
qt = zeros(d,n);
for j=1:n
    [theta, r] = cart2pol(q(1,j),q(2,j));
    R(1,:)=[cos(-theta),-sin(-theta)];
    R(2,:)=[sin(-theta),cos(-theta)];
    k(:,j)=R*k(:,j);
    qtmp = [r,0];
    Fq(:,j)=Flambda(qtmp,a);
    dFk(:,j)=dFlambda(qtmp,k(:,j),a);
    Fqt(:,j) = Fq(:,j) + t*dFk(:,j);
    qt(:,j)=Flambda(Fqt(:,j),1/a);
    [thetat, rt] = cart2pol(qt(1,j),qt(2,j));
    qt(:,j) = [rt*cos(thetat+theta),rt*sin(thetat+theta)];
end
ct = srvf_to_plf(qt,T);


end
