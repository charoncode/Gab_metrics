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


% Plots the geodesic represented by G,T.  G is a 3-dimensional array 
% containing the SRVFs along the geodesic.  G(:,:,k) represents the kth
% SRVF, and all SRVFs have the same change point parameters.
%
% This function plots to the current plot window.  If you want to see 
% all of the functions at the same time, you'll need to set hold on before 
% calling this routine.
%
% This routine can either plot the SRVFs themselves, or it can plot 
% the corresponding curves.  The default is to plot curves; set the argument 
% plot_type to 'q' to draw SRVFs instead.
%
% The color of each curve (or SRVF) along the geodesic is determined by 
% the argument colors, which is a cell array containing at least one color 
% string recognized by the plot() function.  The default is 
% { 'k', 'b', 'c', 'm', 'r' }.  The first curve / SRVF gets colors{1}, 
% the second gets colors{2}, and so on.  Repeats cyclically if there are 
% more curves / SRVFs than colors.
%
% Inputs
%  G : the SRVFs.  G(:,:,k) represents the kth SRVF.
%  T : the change point parameters for the SRVFs in G
%  plot_type : 'f' for functions (the default), or 'q' for SRVFs
%  colors : a cell array containing color strings.  Default is 
%           colors = { 'k', 'b', 'c', 'm', 'r' }.
% --------------------------------------------------------------------------
function plot_geodesic( G, T, plot_type, colors )
  if ( nargin < 3 )
    plot_type = 'f';  % default: plot functions
  end
  if ( nargin < 4 )
    colors = { 'k', 'b', 'c', 'm', 'r' }; % default colors
  end

  [dim nsegs nsteps] = size(G);

  figure();
  hold on

  xl=0; xu=0;
  yl=0; yu=0;
  zl=0; zu=0;

  for i=1:nsteps
    cidx = mod( i-1, length(colors) ) + 1;  % cycle through the colors
    if (i==1)
      %keystr = sprintf('%s;F1;', colors{cidx});
      keystr='b;F1;';
    elseif (i==nsteps )
      %keystr = sprintf('%s;F2;', colors{cidx});
      keystr='r;F2;';
    else
      %keystr = colors{cidx};
      keystr='k';
    end

    if ( plot_type == 'q' )
      Gi = G(:,:,i);
      Gi(1,:) = Gi(1,:) + xu;
      plot_srvf( Gi, T, keystr );

      if ( dim > 1 )
        xl=min(xl, min(Gi(1,:)));
        xu=max(xu, max(Gi(1,:)));
        yl=min(yl, min(Gi(2,:)));
        yu=max(yu, max(Gi(2,:)));
      end
      if ( dim > 2 )
        zl=min(zl, min(Gi(3,:)));
        zu=max(zu, max(Gi(3,:)));
      end
    else
      GFi = srvf_to_plf( G(:,:,i), T );
      GFi(1,:) = GFi(1,:) + xu;
      plot_plf( GFi, T, keystr );

      if ( dim > 1 )
        xl=min(xl, min(GFi(1,:)));
        xu=max(xu, max(GFi(1,:)));
        yl=min(yl, min(GFi(2,:)));
        yu=max(yu, max(GFi(2,:)));
      end
      if ( dim > 2 )
        zl=min(zl, min(GFi(3,:)));
        zu=max(zu, max(GFi(3,:)));
      end
    end
  end

  if ( dim == 1 )
    axis equal
  elseif ( dim == 2 )
    axis([xl xu yl yu], 'equal');
  elseif ( dim == 3 )
  [xl xu yl yu zl zu]
    axis([xl xu yl yu zl zu], 'equal');
  end

  hold off
end
