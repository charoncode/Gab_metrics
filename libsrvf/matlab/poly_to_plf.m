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


% Returns a piecewise-linear function representing the given polygon.
%
% Polygons with dimension 1 are used to represent real-valued functions.  
% For these, this routine produces a PLF which 
%  1) is a reparametrization of the original function P, and 
%  2) always has slope +1 or -1
% Parameter values T for the resulting PLF are calculated automatically, with 
% T(end) equal to the total variation of the function.
% 
% Polygons with dimension 2 or higher are used to represent curves. 
% In this case, the optional argument tv may be used to specify the parameter 
% values corresponding to the vertices in P.  If tv is omitted, then parameter 
% values will be calculated automatically so that, for i=1,2,...,size(P,2), 
% tv(i) = sum of lengths of first (i-1) segments of the polygon.
%
% Inputs:
%  P - the polygon.  An n x k matrix, where n is the dimension of the 
%      ambient space, and k is the number of vertices.
%  tv - parameter values corresponding to the vertices in P.  A 1 x k matrix.
%      Optional, and ignored if size(P,1)==1.
% 
% Outputs:
%  F - the PLF function values.
%  T - the parameter values corresponding to F.
% --------------------------------------------------------------------------
function [F T] = poly_to_plf(P,tv)
  [dim nsamps] = size(P);

  if ( dim > 1 )
    if ( nargin > 1 )
      assert(size(tv,1)==1 && size(tv,2)==size(P,2));
      T = tv;
    else
      T = poly_arclength_params(P);
    end
      
    F = P;
    if ( nargout > 2 )
      error 'S is only useful for 1-D curves';
    end
  else
    S = updown_sequence(P);
    F = [0 cumsum(S)];
    T = [0 cumsum(abs(S))];
  end
end


%!test
%! P=[0 1 1 0; 0 0 1 1];
%! [F T]=poly_to_plf(P);
%! assert(F,P);
%! assert(T,[0 1 2 3]);
%!
%!test
%! P=rand(3,200);
%! tv=linspace(0,1,200);
%! [F T]=poly_to_plf(P,tv);
%! assert(F,P);
%! assert(T,tv);
%!
%!test
%! P=[0 .2 .3 .4 .5 .6 1 .2 0];
%! [F T]=poly_to_plf(P);
%! assert(F,[0 1 0]);
%! assert(T,[0 1 2]);
