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


% Returns a PLF which is equivalent to F,T, but has changepoint parameters 
% Tr.  Tr should be a refinement of T (i.e. T(1)=Tr(1), T(end)=Tr(end), 
% and T is a subset of Tr).
%
% Inputs:
%  F,T - the sample points and changepoint parameters of the PLF
%  Tr - the new changepoint parameter vector
% 
% Outputs:
%  Fr - the sample points corresponding to Tr
% --------------------------------------------------------------------------
function Fr = plf_refine(F,T,Tr)
  Fr = plf_evaluate(F,T,Tr);
end


%!test
%! F=[0 1 0 1 0];
%! T=[0 1/4 1/2 3/4 1];
%! Tr=[0 1/8 1/4 3/8 1/2 5/8 3/4 7/8 1];
%! Fr=plf_refine(F,T,Tr);
%! assert(Fr,[0 1/2 1 1/2 0 1/2 1 1/2 0]);
