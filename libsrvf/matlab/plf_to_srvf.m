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


% Compute the SRVF for a piecewise-linear function.  Since SRVF's are 
% only defined for absolutely-continuous functions, T must be strictly 
% increasing.
%
% Inputs:
%  F - function values
%  T - corresponding parameter values
%
% Outputs:
%  Q - SRVF values.  Q(i) = function value on interval [T(i), T(i+1)]
% --------------------------------------------------------------------------
function Q = plf_to_srvf( F, T )
  assert( min(diff(T)) > 0 );

  [dim, ncp] = size(F);

  if ( dim > 1 )
    V = diff(F,1,2);

    for i=1:dim
      V(i,:) = V(i,:) ./ diff( T );
    end

    Vrmag = sqrt( sqrt( sum( V .* V, 1 ) ) );
    zidx = find( Vrmag < 1e-4 );

    Q = zeros( dim, ncp-1 );
    for i=1:dim
      Q(i,:) = V(i,:) ./ Vrmag;
      Q(i,zidx) = 0;
    end

  else
    m = diff(F) ./ diff(T);
    Q = sign(m) .* sqrt(abs(m));
  end
end


%!test
%! F=[0 1 1 0; 
%!    0 0 1 1];
%! T=[0 1/3 2/3 1];
%! Qexp=[sqrt(3)  0.00000 -sqrt(3);
%!       0.00000  sqrt(3)  0.00000];
%! Q=plf_to_srvf(F,T);
%! assert(Q,Qexp,1e-3);
%!
%!test
%! F=[0 1 0 -1;
%!    0 1 2  1];
%! T=[0 1 2 3];
%! x=2**(-1/4);
%! Qexp=[x -x -x;
%!       x  x -x];
%! Q=plf_to_srvf(F,T);
%! assert(Q,Qexp,1e-3);
%!
%!test
%! F=[0 1 1 0];
%! T=[0 1/3 2/3 1];
%! Qexp=[sqrt(3) 0 -sqrt(3)];
%! Q=plf_to_srvf(F,T);
%! assert(Q,Qexp,1e-5);
