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


% Computes the arclength of the given curve.
%
% For closed curves, F(:,1)==F(:,end) (i.e. the seedpoint is duplicated), 
% so computing the arclength is the same as for open curves.
% 
% Inputs:
%  F,T - the sample points and parameter values of the curve
% 
% Outputs:
%  l - the arclength of the curve
%  S - the arclength function of the curve.  That is, for each i, 
%      S(i) = the arc length of F(:,1) through F(:,i)
% --------------------------------------------------------------------------
function [l S] = plf_arclength( F, T )
  [dim npts] = size( F );

  if ( dim > 1 )
    dF = diff(F,1,2);

    S = [0 cumsum( sqrt(sum( dF .* dF, 1 )), 2 )];
    l = S(end);
  else
    S = [0 cumsum( abs( diff(F) ) )];
    l = S(end);
  end
end


%!test
%! F=[0 1/4 0 1/4 0];
%! T=linspace(0,1,5);
%! assert(plf_arclength(F,T),1,1e-5);
%!
%!test
%! F=[4 -1 -3 5 7 1 0; 3 7 -3 -7 2 1 0; 1 2 3 4 5 6 7];
%! T=0:6;
%! assert(plf_arclength(F,T),42.898,1e-3);

