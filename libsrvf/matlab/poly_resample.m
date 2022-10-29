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


% Given a polygon P in R^n, return a polygon Pr whose vertices are at 
% evenly-spaced intervals along P.
%
% If the optional parameter closed is set to true, then P will be 
% interpreted as a closed polygon.  In that case, the following convention 
% is used regarding the first and last vertices.  If the first vertex of P 
% is equal to the last vertex of P, then the first vertex of Pr will likewise 
% be equal to the last vertex of Pr.  Otherwise, the space between the first 
% and last vertices of Pr will be the same as the spacing between all of the 
% other pairs of adjacent vertices of Pr.
%
% Inputs
%  P - the original polygon: an n x k matrix, where n is the dimension of the 
%      ambient space and k is the number of vertices.
%  nsamps - the number of sample points desired for the result
%  closed - true iff P should be interpreted as a closed polygon
%           Optional, default is false.
%
% Outputs
%  Pr - the resampled polygon
% --------------------------------------------------------------------------
function Pr = poly_resample(P, nsamps, closed)
  % corner cases
  if ( size(P,2) < 2 || nsamps < 2 )
    nsamps=min(nsamps,size(P,2));
    Pr=P(:,1:nsamps);
    return;
  end

  % default: closed=false
  if ( nargin < 3 )
    closed=0;
  end

  % for closed curves, ensure that first vertex = last vertex
  if ( closed )
    if ( norm(P(:,1)-P(:,end)) < 1e-2 )
      % snap closed
      P(:,end)=P(:,1);
      dupseed=1; % seedpoint should be copied
    else
      % add copy of first vertex to the end
      P = [P P(:,1)];
      dupseed=0;  % seedpoint should not be copied
    end
  end

  % compute chord-length parameter values and uniform parameter values
  dP = diff(P,1,2);
  clparams = [0 cumsum(sqrt(sum(dP.*dP,1)))];
  if ( closed && ~dupseed )
    uparams = linspace(0,clparams(end),nsamps+1)(1:nsamps);
  else
    uparams = linspace(0,clparams(end),nsamps);
  end

  % resample
  dim = size(P,1);
  Pr = zeros(dim, nsamps);
  for i=1:dim
    Pr(i,:) = interp1(clparams, P(i,:), uparams);
  end
end


%!demo
%! Po=[0 1 0 -1; 0 1 2 1];
%! Pc=[Po Po(:,1)];
%! Pro=poly_resample(Po,7,false);
%! Prc=poly_resample(Po,8,true);
%! Prcd=poly_resample(Pc,8,true);
%! figure();
%! plot(Po(1,:),Po(2,:),'b-x;Original;',Pro(1,:),Pro(2,:),'r*;Resampled;');
%! title('Open Polygonal Arc, Resampled at 7 Points');
%! figure();
%! plot(Po(1,:),Po(2,:),'b-x;Original;',Prc(1,:),Prc(2,:),'r*;Resampled;');
%! title('Closed Polygon, Resampled at 8 Points, no Duplicate Seedpoint');
%! figure();
%! plot(Pc(1,:),Pc(2,:),'b-x;Original;',Prcd(1,:),Prcd(2,:),'r*;Resampled;');
%! title('Closed Polygon, Resampled at 8 Points, with Duplicate Seedpoint');
