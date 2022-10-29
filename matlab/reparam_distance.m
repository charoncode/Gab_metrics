%% Compute the geodesic distance between two curves using DP or exact algorithm

function d = reparam_distance(c1, T1, c2, T2, method, seed_search, abconst)
    if nargin < 5
        method = 'dp';
        abconst = 1;
        seed_search = 0;
    end
    
    % Checking input arguments
    assert( size(c1,1) == size(c2,1) );
    assert( size(T1,1) == 1 && size(T1,2) == size(c1,2) );
    assert( size(T2,1) == 1 && size(T2,2) == size(c2,2) );
    assert( min(diff(T1)) > 0 );
    assert( min(diff(T2)) > 0 );
    assert( abs(T1(1)-T2(1)) < 1e-12 );
    assert( abs(T1(end)-T2(end)) < 1e-12 );

    [c1n,T1n,c2n,T2n] = curveAlignment(c1,T1,c2,T2,abconst,seed_search,1,method);
    Q1 = SRV_transform(c1n, T1n);
    Q2 = SRV_transform(c2n, T2n);

    d = srvf_abdistance(Q1, T1n, Q2, T2n, abconst);
    

end
