function [G1, GT1, G2, GT2] = optimal_reparam(c1, T1, c2, T2, method, abconst)

    if strcmp(method, 'klr') % Exact matching following the paper by
                             % Klassen, Lahiri and Robinson 

        % Compute SRVFs
        Q1 = SRV_transform(c1, T1);
        Q2 = SRV_transform(c2, T2);

        % Find optimal matching
        [G1, GT1, G2, GT2] = ...
            srvf_klr_optimal_reparam(Q1, T1, Q2, T2, abconst);     
        
    elseif strcmp(method, 'dp') % Dynamic programming
       [G1, GT1] = dynamic_programming(c1, T1, c2, T2, abconst);
       N = size(G1,2);
       G2 = linspace(0,1,N);
       GT2 = linspace(0,1,N);
    
    else
        error('Unknown method.');
    end