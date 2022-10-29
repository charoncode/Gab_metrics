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


% Computes the L^2 distance between two SRVF's.
% The SRVFs must be defined on the same interval.
%
% Inputs
%  Q1,T1 : the first SRVF
%  Q2,T2 : the second SRVF
%
% Outputs
%  d : the L^2 distance between the SRVFs
% --------------------------------------------------------------------------
function d = srvf_l2distance( Q1, T1, Q2, T2 )
  [Qd Td] = srvf_linear_combination( Q1, T1, Q2, T2, 1, -1 );
  d = srvf_l2norm( Qd, Td );
end


%!test
%! Q1=[1 -1];
%! T1=[0 1/4 1];
%! Q2=[-1 1 2];
%! T2=[0 1/2 3/4 1];
%! assert(srvf_l2distance(Q1,T1,Q2,T2),sqrt(4.25),1e-4);
