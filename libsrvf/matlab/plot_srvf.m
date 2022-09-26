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


% Plots the given SRVF.
%
% Inputs
%  Q,T : the sample points and changepoint parameter values of the given SRVF
%  fmt : a format string recognized by the plot() function
% --------------------------------------------------------------------------
function plot_srvf( Q, T, fmt )
  if ( nargin < 3 )
    fmt = 'b';
  end

  [dim nsegs] = size(Q);
  
  if ( dim == 1 )
    hold on;
    for i=1:length(Q)
      H=plot([T(i) T(i+1)], [Q(i), Q(i)], fmt );
      set(H,'linewidth',2);
    end
  elseif ( dim == 2 )
    H=plot(Q(1,:),Q(2,:),fmt);
  elseif ( dim == 3 )
    H=plot3(Q(1,:),Q(2,:),Q(3,:),fmt);
  else
    error( 'SRVFs must be dimension 3 or lower.' );
  end
end
