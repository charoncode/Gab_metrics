function [gam, Tr] = dynamic_programming(c1i,T1i,c2i,T2i,a)

N  = size(c1i,2);
T1 = linspace(0,1,N);
T2 = linspace(0,1,N);
Tr = T1;

[T, G1] = plf_constant_speed_reparam(c1i,T1i);
c1 = plf_evaluate(c1i,G1,T1);

[T, G2] = plf_constant_speed_reparam(c2i,T2i);
c2 = plf_evaluate(c2i,G2,T2);

Q1 = SRV_transform(c1, T1);
Q2 = SRV_transform(c2, T2);

gamtm = DynamicProgrammingQ(Q1,Q2,a,0,0);
gam = zeros(1,size(gamtm,2)+1);
gam(2:end) = gamtm;
gam = (gam-gam(1))/(gam(end)-gam(1));