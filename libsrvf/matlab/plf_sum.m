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


% Returns the pointwise sum of two PLFs.
%
% Inputs:
%  F1,T1 - the sample points and parameter values of the first PLF
%  F2,T2 - the sample points and parameter values of the second PLF
%
% Outputs:
%  F,T - the sample points and parameter values of the sum PLF
% --------------------------------------------------------------------------
function [F T] = plf_sum( F1, T1, F2, T2 )
  T = unique( [T1 T2] );
  F = plf_evaluate( F1, T1, T ) + plf_evaluate( F2, T2, T );
end


%!test
%! F1 = [0 1 2 3 4];
%! T1 = [0 1/4 1/2 3/4 1];
%! F2 = [5 6 7 8 9];
%! T2 = [0 1/8 3/8 5/8 1];
%! Texp = [0 1/8 1/4 3/8 1/2 5/8 3/4 1];
%! Fexp = [5.0 6.5 7.5 8.5 9.5 10.5 11.3333 13.0];
%! [F T] = plf_sum(F1,T1,F2,T2);
%! assert( F, Fexp, 1e-4 );
%! tv = linspace(0,1,100);
%! F1tv = interp1(T1,F1,tv);
%! F2tv = interp1(T2,F2,tv);
%! Ftv = interp1(T,F,tv);
%! assert( Ftv, F1tv+F2tv, 1e-4 );
