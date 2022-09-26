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


% Computes the L^2 inner product of two SRVFs.  The SRVFs must have the 
% same dimension and must be defined on the same interval.
%
% Inputs
%  Q1,T1 : the first SRVF
%  Q2,T2 : the second SRVF
%
% Outputs
%  ip : the L^2 inner product of the two SRVFs
% ---------------------------------------------------------
function ip = srvf_l2product( Q1, T1, Q2, T2 )
  assert( size(Q1,1) == size(Q2,1) );
  assert( size(T1,1) == 1 );
  assert( size(T2,1) == 1 );
  assert( size(Q1,2) == size(T1,2)-1 );
  assert( size(Q2,2) == size(T2,2)-1 );

  Tr = unique( [T1 T2] );

  Q1r = srvf_refine( Q1, T1, Tr );
  Q2r = srvf_refine( Q2, T2, Tr );

  ip = sum( diff(Tr) .* sum(Q1r.* Q2r,1) );
end


%!test
%! Q1=[0 0 0 0];
%! T1=[0 1/8 1/2 3/4 1];
%! Q2=[1];
%! T2=[0 1];
%! assert( srvf_l2product(Q1,T1,Q2,T2), 0, 1e-6 );
%!
%!test
%! Q1=[1 -1 1 -1 1];
%! T1=linspace(0,1,6);
%! Q2=[1];
%! T2=[0 1];
%! assert( srvf_l2product(Q1,T1,Q2,T2), 0.2, 1e-6 );
