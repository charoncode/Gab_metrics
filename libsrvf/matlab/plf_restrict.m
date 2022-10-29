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


% Returns a new PLF representing the restriction of F,T to [a,b].
% [a,b] should be contained in [T(1),T(end)] (the original domain of the PLF).
%
% Inputs:
%  F,T - the sample points and parameter values of the PLF
%  a,b - the endpoints of the new domain interval [a,b]
%
% Outputs:
%  Fr, Tr - the sample points and parameter values of the restricted PLF
% --------------------------------------------------------------------------
function [Fr Tr] = plf_restrict(F,T,a,b)
  Tr = unique([T a b]);
  Fr = plf_refine(F,T,Tr);
  aidx = find(Tr==a);
  bidx = find(Tr==b);
  Tr = Tr(aidx:bidx);
  Fr = Fr(:,aidx:bidx);
end


%!test
%! F=[0 1 0 1 0; 1 2 1 2 1];
%! T=[0 1/4 1/2 3/4 1];
%! a=1/8;
%! b=7/8;
%! Frx=[1/2 1 0 1 1/2; 3/2 2 1 2 3/2];
%! Trx=[1/8 1/4 1/2 3/4 7/8];
%! [Fr Tr]=plf_restrict(F,T,a,b);
%! assert(Tr,Trx);
%! assert(Fr,Frx);
