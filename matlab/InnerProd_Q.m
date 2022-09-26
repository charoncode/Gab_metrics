function val = InnerProd_Q(q1,T1,q2,T2)

% T1(size(T1,2)) = (1+T1(size(T1,2)-1))/2;
% T2(size(T2,2)) = (1+T2(size(T2,2)-1))/2;
% T1(size(T1,2)+1) = 1;
% T2(size(T2,2)+1) = 1;
T = unique( [T1 T2] );
q1r=plf_refine( q1, T1, T );
q2r=plf_refine( q2, T2, T );
%[~,T] = size(q1r);

val = trapz(T,sum(q1r.*q2r));
