function k=deriv_srvf(c,h,T)

deriv_c = diff(c,1,2)./diff(T);
deriv_h = diff(h,1,2)./diff(T);
proj_tang = srvf_l2product(deriv_c,T,deriv_h,T)*deriv_c/srvf_l2norm(deriv_c,T);
proj_orth = deriv_h - proj_tang;
norm_c = srvf_l2norm(deriv_c,T);
k = 1/(sqrt(norm_c)) * (proj_orth+1/2*proj_tang);

end
