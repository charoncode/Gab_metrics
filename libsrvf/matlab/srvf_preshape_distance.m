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


% Computes the great circle distance between two SRVFs.
%
% Inputs:
%  Q1,T1 : the sample points and parameter values of the first SRVF
%  Q2,T2 : the sample points and parameter values of the second SRVF
%
% Outputs:
%  d : the great circle distance from Q1 to Q2 on the sphere
% --------------------------------------------------------------------------
function d = srvf_preshape_distance( Q1, T1, Q2, T2 )
  ip = srvf_l2product( Q1, T1, Q2, T2 );
  n1 = srvf_l2norm( Q1, T1 );
  n2 = srvf_l2norm( Q2, T2 );

  % Preshape distance isn't defined if SRVF's have different norms.
  assert( abs(n1-n2) < max(0.01*max(n1,n2), 1e-3) );
  
  if ( n1*n2 > 1e-4 )
    x = ip / n1 / n2;
    if ( x < -1.001 || x > 1.001 )
      error('srvf_preshape_distance: x=%f out of bounds\n', x);
    elseif ( x < -1 )
      x=-1;
    elseif ( x > 1 )
      x=1;
    end
    d=acos(x);
  else
    d = 0;
  end

end
