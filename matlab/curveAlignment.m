function [c1new,T1new,c2new,T2new] = curveAlignment(c1i,T1i,c2i,T2i,a,seed_search_flag,scaling_flag,method)

if nargin == 5
    seed_search_flag = 0;
    scaling_flag = 1;
end
if nargin == 6
    scaling_flag = 1;
end


N  = size(c1i,2);
T1 = linspace(0,1,N);
T2 = linspace(0,1,N);

[T, G1] = plf_constant_speed_reparam(c1i,T1i);
c1 = plf_evaluate(c1i,G1,T1);

[T, G2] = plf_constant_speed_reparam(c2i,T2i);
c2 = plf_evaluate(c2i,G2,T2);


c1 = preprocess_curve(c1,scaling_flag);
c2 = preprocess_curve(c2,scaling_flag);


% Extra step for closed curves (find the starting point)

if seed_search_flag == 1
    % perturb curves if necessary to avoid an error
    if c1(:,N) == c1(:,1)
        c1(:,N) = .99*c1(:,N)+.01*c1(:,N-1);
    end

    if c2(:,N) == c2(:,1)
        c2(:,N) = .99*c2(:,N)+.01*c2(:,N-1);
    end
    
% Find best seed
dists = zeros(1,N);

for j = 1:N
    c2n = ShiftF(c2,j);
    c2n=preprocess_curve(c2n, scaling_flag);
    Q1 = SRV_transform(c1, T1);
    Q2n = SRV_transform(c2n, T2);
    R = srvf_optimal_rotation( Q1, T1, Q2n, T2 );
    c2n = R*c2n;
    Q2n = SRV_transform(c2n, T2);
    Q1n = SRV_transform(c1, T1);
    Q2n = SRV_transform(c2n, T2);
    dists(j) = srvf_abdistance( Q1n, T1, Q2n, T2, 1 );
end

[~,sorted_inds]=sort(dists);


c2 = ShiftF(c2,sorted_inds(1));
c2 = preprocess_curve(c2,scaling_flag);
end

% Send the curves to SRV transform space
Q1 = SRV_transform(c1, T1);
Q2 = SRV_transform(c2, T2);

% Register transformed curves over rotations
R = srvf_optimal_rotation( Q1, T1, Q2, T2 );
c2 = R*c2;
Q2 = SRV_transform(c2, T2);

% Find optimal reparameterization
[G1, GT1, G2, GT2] = optimal_reparam(c1, T1, c2, T2, method, a);

% Apply optimal matching
c1 = Group_Action_by_Gamma_Coord(c1,G1);
c2 = Group_Action_by_Gamma_Coord(c2,G2);

c1new = c1;
T1new = GT1;
c2new = c2;
T2new = GT2;
