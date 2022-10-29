% libsrvf
% =======
%
% A shape analysis library using the square root velocity framework.
% 
% Copyright (C) 2012   FSU Statistical Shape Analysis and Modeling Group
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>
% --------------------------------------------------------------------


% Computes the covariance of a collection of points on the unit sphere in L^2.
% We first lift the points into the tangent space at Qm, and then we compute 
% the covariance of the resulting collection of points in L^2.
%
% Inputs:
%  Qs : the collection of SRVFs.  Must all have the same domain, and 
%       must already be aligned to Qm.
%  Ts : the corresponding parameter values for Qs
%  Qm : the Karcher mean of the functions in Qs.  Must have same domain as 
%       the SRVFs in Qs.
%  Ts : the corresponding parameter values for Qm
% --------------------------------------------------------------------------
function [K, Mu] = srvf_covariance(Qs, Ts, Qm, Tm)
  uv = linspace(Tm(1), Tm(end), length(Tm));

  % TODO: to get a bunch of vectors of the same length, we're evaluating 
  % all of the SRVF's at uniformly-spaced parameters.  Almost certainly 
  % not the best way to do this.
  Dat = zeros(length(Qs), size(Qm,1)*length(uv));
  for i=1:length(Qs)
    [V TV] = sphere_shooting_vector(Qm, Tm, Qs{i}, Ts{i});
    Vsamps = srvf_evaluate(V, TV, uv);
    Dat(i,:) = Vsamps'(:);
  end

  K = cov(Dat);
  Mu = mean(Dat);
end


%!demo
%! debug_on_error(1);
%! load ../tests/data/rna1.mat
%! load ../tests/data/rna2.mat
%! load ../tests/data/rna3.mat
%! load ../tests/data/rna4.mat
%! load ../tests/data/rna5.mat
%! nfuncs=5;
%! colors = {'b', 'g', 'c', 'm', 'r'};
%! 
%! [Fs{1}, Ts{1}] = poly_to_plf(X1);
%! [Fs{2}, Ts{2}] = poly_to_plf(X2);
%! [Fs{3}, Ts{3}] = poly_to_plf(X3);
%! [Fs{4}, Ts{4}] = poly_to_plf(X4);
%! [Fs{5}, Ts{5}] = poly_to_plf(X5);
%! figure();
%! hold on;
%! for i=1:nfuncs
%!   Ts{i} = Ts{i} / Ts{i}(end);
%!   Qs{i} = plf_to_srvf(Fs{i}, Ts{i});
%!   nrm = srvf_l2norm(Qs{i}, Ts{i});
%!   Qs{i} = Qs{i} / nrm;
%!   Fs{i} = Fs{i} / (nrm*nrm);
%!   plot3(Fs{i}(1,:), Fs{i}(2,:), Fs{i}(3,:), colors{i});
%! end
%! [Qm, Tm] = srvf_karcher_mean(Qs, Ts, 1e-3, 30, 0.3);
%! Fm = srvf_to_plf(Qm, Tm);
%! figure();
%! hold on;
%! plot3(Fm(1,:), Fm(2,:), Fm(3,:), 'k;Fm;');
%! for i=1:nfuncs
%!   R = srvf_optimal_rotation(Qm, Tm, Qs{i}, Ts{i});
%!   Qs{i} = R*Qs{i};
%!   [G T] = srvf_optimal_reparam(Qm, Tm, Qs{i}, Ts{i});
%!   [Qs{i}, Ts{i}] = srvf_gamma_action(Qs{i}, Ts{i}, G, T);
%!   R = srvf_optimal_rotation(Qm, Tm, Qs{i}, Ts{i});
%!   Qs{i} = R*Qs{i};
%!   Fsr{i} = srvf_to_plf(Qs{i}, Ts{i});
%!   plot3(Fsr{i}(1,:), Fsr{i}(2,:), Fsr{i}(3,:), colors{i});
%! end
%! [K, Mu] = srvf_covariance(Qs, Ts, Qm, Tm)
