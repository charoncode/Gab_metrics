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


% Returns the PLF corrresponding the the given SRVF.  The output PLF 
% will have the same change point vector as the input SRVF, so no change 
% point vector is returned.
%
% Inputs
%  Q,T : the SRVF
%
% Outputs
%  F : the PLF function values
% --------------------------------------------------------------------------
function F = srvf_to_plf( Q, T )
  [dim nsamps] = size(Q);

  if ( dim > 1 )
    normQ = sqrt( sum( Q.*Q, 1 ) );
    V = Q .* repmat( normQ, dim, 1 ) .* repmat( diff(T), dim, 1 );
    F = [zeros(dim,1) cumsum( V, 2 )];
  else
    v = Q .* abs(Q);
    F = [0 cumsum(v .* diff(T))];
  end
end


%!assert(srvf_to_plf([0],[0 1]),[0 0],1e-6);
%!assert(srvf_to_plf([0;0],[0 1]),[0 0;0 0],1e-6);
%!assert(srvf_to_plf([1],[0 1]),[0 1],1e-6);
%!assert(srvf_to_plf([1;1],[0 1]),[0 sqrt(2);0 sqrt(2)],1e-6);
%!assert(srvf_to_plf([1 -1],[0 1/2 1]),[0 1/2 0],1e-6);
%!assert(srvf_to_plf([1 -1;0 1],[0 1/2 1]),[0 1/2 -.20711;0 0 sqrt(2)/2],1e-4);
