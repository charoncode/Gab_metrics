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


% Returns the arclength parameter values for a polygon P.
%
% Inputs:
%  P - the polygon.  An n x k matrix, where n is the dimension of the 
%      ambient space, and k is the number of vertices.
%
% Outputs:
%  clp - clp(i) = sum of lengths of the first (i-1) segments of P
% --------------------------------------------------------------------------
function clp = poly_arclength_params(P)
  dP = diff(P,1,2);
  clp = [0 cumsum(sqrt(sum(dP.*dP,1)))];
end


%!test
%! P=[0 1 1 0; 0 0 1 1];
%! clp=poly_arclength_params(P);
%! assert(clp,[0 1 2 3]);
