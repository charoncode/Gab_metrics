function d = srvf_abdistance( Q1, T1, Q2, T2, a )
  T = unique( [T1 T2] );

  Q1r = srvf_refine( Q1, T1, T );
  Q2r = srvf_refine( Q2, T2, T );

  theta = real(acos(dot(Q1r,Q2r)./(sqrt(dot(Q1r,Q1r)).*sqrt(dot(Q2r,Q2r)))));
  W = (sqrt(dot(Q1r,Q1r).*diff(T)).*sqrt(dot(Q2r,Q2r).*diff(T))).*((a*theta<pi/2).*cos(a*theta));
  TF = isnan(W);
  W(TF) = 0;
  prod = sum(W);
  len = (sum(sum(Q1r.*Q1r,1).*diff(T) + sum(Q2r.*Q2r,1).*diff(T)));
  d = sqrt(len - 2*prod);
end

