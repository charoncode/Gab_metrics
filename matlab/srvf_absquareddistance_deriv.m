function deriv = srvf_absquareddistance_deriv(Q1,T1,Q2,T2,a)

N = size(Q1,2);

T = unique( [T1 T2] );

Q1r = srvf_refine( Q1, T1, T );
Q2r = srvf_refine( Q2, T2, T );

theta = real(acos(dot(Q1r,Q2r)./(sqrt(dot(Q1r,Q1r)).*sqrt(dot(Q2r,Q2r)))));
W = (sqrt(dot(Q1r,Q1r).*diff(T)).*sqrt(dot(Q2r,Q2r).*diff(T))).*((a*theta<pi/2).*sin(a*theta).*theta);
TF = isnan(W);
W(TF) = 0;
prod = sum(W);
deriv = 2*prod;