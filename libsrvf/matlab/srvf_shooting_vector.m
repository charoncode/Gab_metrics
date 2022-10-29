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


% Returns the shooting vector from Q1 to Q2.
%
% Given two points Q1 and Q2 on the same sphere in L^2, returns a vector V 
% in the tangent space at Q1 which is tangent to a geodesic from Q1 to Q2 
% and has length equal to the geodesic distance from Q1 to Q2.
%
% Inputs:
%  Q1,T1 : the sample points and parameter values of the first point
%  Q2,T2 : the sample points and parameter values of the second point
% 
% Outputs:
%  V,Tv : the sample points and parameter values of the shooting vector
% --------------------------------------------------------------------------
function [V TV] = srvf_shooting_vector(Q1, T1, Q2, T2)
  Q1nrm = srvf_l2norm(Q1, T1);
  proj_nrm = srvf_l2product(Q1, T1, Q2, T2) / (Q1nrm*Q1nrm);
  d = srvf_preshape_distance(Q1, T1, Q2, T2);
  [V TV] = srvf_linear_combination(Q2, T2, Q1, T1, 1, -proj_nrm);
  Vnrm = srvf_l2norm(V, TV);
  if (Vnrm > 1e-6)
    V = (d/Vnrm)*V;
  else
    V = zeros(size(Q1,1), length(TV)-1);
  end
end


%!test
%! Q1 = [1];
%! T1 = [0 1];
%! Q2 = [1 -1];
%! T2 = [0 0.5 1];
%! Vexp = (pi/2)*[1 -1];
%! TVexp = [0 0.5 1];
%! [V TV] = srvf_shooting_vector(Q1, T1, Q2, T2);
%! assert(V,Vexp);
%! assert(TV,TVexp);

