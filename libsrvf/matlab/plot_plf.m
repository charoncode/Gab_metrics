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


% Plots the given PLF.
%
% Inputs
%  F,T : the sample points and parameter values of the given PLF
%  fmt : a format string recognized by the plot() function
% --------------------------------------------------------------------------
function plot_plf( F, T, color )
  if ( nargin < 3 )
    color = 'b';
  end

  [dim nsegs] = size(F);
  
  if ( dim == 1 )
    plot(T,F,color);
  elseif ( dim == 2 )
    plot(F(1,:),F(2,:),color);
  elseif ( dim == 3 )
    plot3(F(1,:),F(2,:),F(3,:),color);
  else
    error( 'PLFs must be dimension 3 or lower.' );
  end
end
