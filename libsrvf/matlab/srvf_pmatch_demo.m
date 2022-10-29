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


% A simple demo of the libsrvf partial matching functions.
%
% Inputs:
%  X1 - sample points for the first curve
%  X2 - sample points for the second curve
%  W - grid width: the number of break points for X1 (default = 30)
%  H - grid height: the number of break points for X2 (default = 30)
%  do_rots - nonzero to optimize over rotations
% --------------------------------------------------------------------------
function srvf_pmatch_demo(X1,X2,W,H,do_rots)
  if (nargin < 3) W=30; end
  if (nargin < 4) H=30; end
  if (nargin < 5) do_rots=1; end

  % Assign uniformly-spaced parameter values to sample points
  T1 = linspace(0,1,size(X1,2));
  T2 = linspace(0,1,size(X2,2));

  % Scale curves down by a common factor
  L1 = plf_arclength(X1,T1);
  L2 = plf_arclength(X2,T2);
  L = min(L1,L2);
  X1 = X1 / L;
  X2 = X2 / L;

  % Find the Pareto set
  P = srvf_pmatch_find_matches(X1,T1,X2,T2,W,H,do_rots);

  % Select a single match from the Pareto set
  sel_idx = srvf_pmatch_select_single(T1,T2,P,0.01)
  ai = P(sel_idx,1); bi = P(sel_idx,2);
  ci = P(sel_idx,3); di = P(sel_idx,4);
  X1sel = X1(:,ai:bi);
  X2sel = X2(:,ci:di);
  T1sel = linspace(0,1,size(X1sel,2));
  T2sel = linspace(0,1,size(X2sel,2));

  % Translate subcurves to the origin, rotationally align them, 
  % and display the match
  X1sel = plf_trans_to_origin(X1sel,T1sel);
  X2sel = plf_trans_to_origin(X2sel,T2sel);
  Q1sel = plf_to_srvf(X1sel,T1sel);
  Q2sel = plf_to_srvf(X2sel,T2sel);
  R = srvf_optimal_rotation(Q1sel,T1sel,Q2sel,T2sel);
  X2sel = R*X2sel;

  figure();
  hold on;
  plot_plf(X1,T1,'b');
  plot_plf(X2,T2,'r');
  title 'Original curves';
  axis equal;

  figure();
  hold on;
  plot_plf(X1sel,T1sel,'b');
  plot_plf(X2sel,T2sel,'r');
  axis equal;
  title 'Selected match';
end
