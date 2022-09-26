function Q = SRV_transform( F, T )
  assert( min(diff(T)) > 0 );

  [dim, ncp] = size(F);

  if ( dim > 1 )
    V = diff(F,1,2);

    for i=1:dim
      V(i,:) = V(i,:) ./ diff( T );
    end

    Vrmag = sqrt( sqrt( sum( V .* V, 1 ) ) );
    zidx = find( Vrmag < 1e-4 );

    Q = zeros( dim, ncp-1 );
    for i=1:dim
      Q(i,:) = V(i,:) ./ Vrmag;
      Q(i,zidx) = 0;
    end

  else
    m = diff(F) ./ diff(T);
    Q = sign(m) .* sqrt(abs(m));
  end
end