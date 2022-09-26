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


% Returns a PLF which is an inverse of the given PLF.  The given PLF 
% must be 1-dimensional and must be monotone (i.e. non-decreasing or 
% non-increasing).  If F is not strictly monotone, then the corresponding 
% function is not one-to-one, and the result will be a right-inverse of 
% the input function.
%
% Inputs:
%  F,T - the sample points and parameter values of the PLF
%
% Outputs:
%  Fi, Ti - the sample points and parameter values of the inverse PLF
% --------------------------------------------------------------------------
function [Fi Ti] = plf_inverse( F, T )
  assert( size(F,1) == 1 ); % 1-D functions only
  assert( ismonotone(F,0) ); % non-increasing or non-decreasing
  assert( min(diff(T)) >= 0 ); % non-decreasing (function may have jumps)

  %if ( nargin < 3 )
  %  min_diff = 1e-8;
  %end

  %Ti = make_monotone( F, min_diff );
  %Fi = make_monotone( T, min_diff );

  %Ti = Ti / sum(diff(Ti));

  Fi = T;
  Ti = F;
end


%function v = make_monotone( v, min_diff )
%  d = guess_direction( v );
%  assert( abs(d) > 0.5 );
%  
%  for i=2:length( v )
%    if ( d > 0 ) 
%      v(i) = max( v(i), v(i-1) + min_diff );
%    else
%      v(i) = min( v(i), v(i-1) - min_diff );
%    end
%  end
%end
%
%function d = guess_direction( v )
%  d = sign( sum( v ) );
%end
