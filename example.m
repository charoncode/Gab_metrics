addpath ./'dynamic programming'/
addpath(genpath('libsrvf\'))
addpath matlab\
addpath metric_estimation\
addpath data\

load prot


%% Choose the different constants and settings
Tr = linspace(0,1,100);


% Choose the two curves (last digit between 1 and 1300)
c1i=prot{1};
c2i=prot{2};
T1=linspace(0,1,size(c1i,2));
T2=linspace(0,1,size(c2i,2));

% Refine the number of curves
c1 = plf_refine(c1i,T1,Tr);
c2 = plf_refine(c2i,T2,Tr);


% Choose the constants for the metric (a=1, b=1/2 for classic srv metric)
a = 2;
b = 1/2;

% Choose the method to compute optimal reparametrisation
 % 'dp' for dynamic programming
 % 'klr' for exact algorithm (can take a long time, use plf_refine to
 % downsample the number of points)
method = 'klr';

% Do you want to compute geodesic ?
compGeod = 'y'; % 'y' for yes, 'n' for no
Nsteps = 5; % number of steps for geodesic

% Do you want to plot the results ?
plotRes = 'y'; 



%% Compute optimal reparametrization

% First, pre-alignment step for rotation and first seed :

[T, G1] = plf_constant_speed_reparam(c1,Tr);
c1 = plf_evaluate(c1,G1,T);

[T, G2] = plf_constant_speed_reparam(c2,Tr);
c2 = plf_evaluate(c2,G2,T);


c1 = preprocess_curve(c1,1);
c2 = preprocess_curve(c2,1);


Q1 = SRV_transform(c1, T);
Q2 = SRV_transform(c2, T);

% Register transformed curves over rotations
R = srvf_optimal_rotation( Q1, T, Q2, T );
c2=R*c2;


[dim, nsegs] = size(c1);

acons = a/(2*b);

[G1,T11,G2,T22] = optimal_reparam(c1,T,c2,T,method,acons);

if strcmp(plotRes, 'y')
    figure
    plot(G1,G2)
end


c1n = plf_evaluate(c1,Tr,G1);
c2n = plf_evaluate(c2,Tr,G2);
[dim, nsegs] = size(c1n);


%% Compute geodesics

if strcmp(compGeod, 'y')
Tt = linspace(0,1,nsegs);
timeGeod = linspace(0,1,Nsteps);
ct = zeros( dim, nsegs, Nsteps );

for i=1:Nsteps
    t = timeGeod(i);
    ct(:,:,i)=geodesic( c1n, Tt, c2n, Tt, t, acons, 3 );

    if strcmp(plotRes, 'y')
    figure
    plot3(ct(1,:,i),ct(2,:,i),ct(3,:,i))
    end
end
end