function [Ft T] = geodesic_deriv( F1, T1, F2, T2, t, abcons )
  assert( size(F1,1) == size(F2,1) );
  assert( min(diff(T1)) > 0 );
  assert( min(diff(T2)) > 0 );
  assert( abs(T1(1)-T2(1)) < 1e-4 );
  assert( abs(T1(end)-T2(end)) < 1e-4 );

  T = unique( [T1 T2] );
  Q1=plf_to_srvf(F1,T1);
  Q2=plf_to_srvf(F2,T2);

  Q1r = srvf_refine( Q1, T1, T );
  Q2r = srvf_refine( Q2, T2, T );
  [dim nsegs] = size(Q1r);

  P = zeros( dim, nsegs );

  ip = srvf_l2product( Q1r, T, Q2r, T );
  theta = acos( ip );

  for j=1:nsegs
    P(:,j) = geodP_deriv(Q1r(:,j), Q2r(:,j), t, abcons);
    %P(:,:,i)=(sin(theta*(1-Ptv(i)))*Q1r + sin(theta*Ptv(i))*Q2r) ./ sin(theta);
  end

  Ft = srvf_to_plf(P,T);
  


