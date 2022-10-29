% LibSRVF - a shape analysis library using the square root velocity framework.
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


% Checks vector v for monotonicity, and returns the direction (increasing, 
% decreasing, constant, or none).
%
% Meaning of parameter tol:
% - if tol==0, check for non-increasing or non-decreasing sequence (default).
% - if tol>0, allow backward steps of size <= tol
% - if tol<0, require forward steps of size >= tol
%
% Inputs
%  v:   vector to check for monotonicity
%  tol: see above
% 
% Outputs
%  b: a bitfield indicating monotonicity.  Can be tested as follows:
%   bitand(b,1)==true  -->  v is increasing (within tolerance)
%   bitand(b,2)==true  -->  v is decreasing (within tolerance)
%   bitand(b,3)==true  -->  v is both increasing and decreasing
%                           (i.e. v is constant, within tolerance).
% --------------------------------------------------------------------------
function b = ismonotone( v, tol )
  if ( nargin < 2 )
    tol = 0;
  end

  b = 0;
  dv = diff(v);
  if ( min(dv) >= -tol ) b = bitor( b, 1 ); end
  if ( max(dv) <= tol ) b = bitor( b, 2 ); end
end


%!test assert(ismonotone(linspace(0,1,20)),1);
%!test assert(ismonotone(linspace(1,0,20)),2);
%!test assert(ismonotone(zeros(1,100)),3);
%!test
%! v=[0 -0.01 0 0 0.01 0.25 1];
%! assert(ismonotone(v,0.011),1);
%! assert(ismonotone(-v,0.011),2);
