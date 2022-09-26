function p = FaTransform(c,T,a)

% Transforms a plane curve c into another plane curve p according to the
% a-transform. 

[d,n]=size(c);

v=curveDeriv(c).*gradient(T,1/n);

% [theta,rho]=curve2Polar(v);
% psi=a*theta;
% nu=sqrt(rho);
% 
% [x(:),y(:)]=pol2cart(psi,nu);

p=zeros(d,n);

% for j=1:n
%     p(1,j)=x(j);
%     p(2,j)=y(j);
% end

p = v/sqrt(sqrt(InnerProd_Q(v,T,v,T)));

