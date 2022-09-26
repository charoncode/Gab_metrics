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


% Computes several steps along the minimum-length great circle path between 
% the given SRVFs Q1,T1 and Q2,T2.  The SRVFs must have the same dimension 
% and the same L^2 norm, and must be defined on the same interval.
%
% If you want a geodesic in the shape space (i.e. between the orbits of 
% Q1 and Q2), then you should first optimally rotate and reparametrize Q2 
% using srvf_optimal_rotation() and srvf_optimal_reparam().
%
% Inputs:
%  Q1,T1  : the SRVF at the beginning of the geodesic
%  Q2,T2  : the SRVF at the end of the geodesic
%  Nsteps : the number of points (including the end points) along 
%           the geodesic to return
%
% Outputs
%  P : a 3-dimensional array containing the points on the geodesic.
%      P(:,:,k) is an SRVF representing the kth point.  All of these 
%      SRVFs will have the same change point parameters.
%  T : a 1-row matrix containing the change point parameters for the SRVFs in P.
% --------------------------------------------------------------------------
function [P T] = srvf_geodesic( Q1, T1, Q2, T2, Nsteps )
  assert( size(Q1,1) == size(Q2,1) );
  assert( size(T1,1) == 1 && size(T1,2) == size(Q1,2)+1 );
  assert( size(T2,1) == 1 && size(T2,2) == size(Q2,2)+1 );
  assert( min(diff(T1)) > 0 );
  assert( min(diff(T2)) > 0 );
  assert( abs(T1(1)-T2(1)) < 1e-4 );
  assert( abs(T1(end)-T2(end)) < 1e-4 );

  T = unique( [T1 T2] );
  Q1r = srvf_refine( Q1, T1, T );
  Q2r = srvf_refine( Q2, T2, T );
  [dim nsegs] = size(Q1r);

  P = zeros( dim, nsegs, Nsteps );
  Ptv = linspace( 0, 1, Nsteps );

  ip = srvf_l2product( Q1r, T, Q2r, T );
  theta = acos( ip );

  for i=1:Nsteps
    P(:,:,i)=(sin(theta*(1-Ptv(i)))*Q1r + sin(theta*Ptv(i))*Q2r) ./ sin(theta);
  end
end


%!demo
%! load ../tests/data/horse-1.csv
%! load ../tests/data/horse-2.csv
%! X1 = horse_1;
%! X2 = horse_2;
%! [F1 T1]=poly_to_plf(X1);
%! [F2 T2]=poly_to_plf(X2);
%! Q1=plf_to_srvf(F1,T1);
%! Q2=plf_to_srvf(F2,T2);
%! [G GT]=srvf_optimal_matching(Q1,T1,Q2,T2);
%! [F2r T2r]=plf_compose(F2,T2,G,GT);
%! Q2r=plf_to_srvf(F2r,T2r);
%! [P T] = srvf_geodesic(Q1,T1,Q2r,T2r,5);
%! plot_geodesic(P,T);
%! title('The curves corresponding to the geodesic');

%!#demo
%! load demos/rna1.mat
%! load demos/rna2.mat
%! [F1 T1]=poly_to_plf(X1);
%! [F2 T2]=poly_to_plf(X2);
%! Q1=plf_to_srvf(F1,T1);
%! Q2=plf_to_srvf(F2,T2);
%! [G GT]=srvf_optimal_matching(Q1,T1,Q2,T2);
%! [F2r T2r]=plf_compose(F2,T2,G,GT);
%! Q2r=plf_to_srvf(F2r,T2r);
%! [P T] = srvf_geodesic(Q1,T1,Q2r,T2r,5);
%! plot_geodesic(P,T);
%! title('The curves corresponding to the geodesic');
%! plot_geodesic(P,T,'q');
%! title('The SRVFs on the geodesic');
