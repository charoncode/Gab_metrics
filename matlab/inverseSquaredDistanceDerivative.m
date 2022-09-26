function deriv = inverseSquaredDistanceDerivative(Q1,T1,Q2,T2,a)

N = size(c1,2);


numerator = -srvf_absquareddistance_deriv(Q1,T1,Q2,T2,a);

denominator = srvf_abdistance( Q1, T1, Q2, T2, a )^2;

deriv = numerator/denominator;