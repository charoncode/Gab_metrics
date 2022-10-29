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


% Translates the sample points of the given PLF so that they have zero mean.
%
% Inputs:
%  F,T - the sample points and parameter values of the given PLF
%
% Outputs:
%  Fc - the centered sample points
% --------------------------------------------------------------------------
function Fc=plf_trans_to_origin(F,T);
  [dim npts]=size(F);
  ctr=eye(npts) - ones(npts) * (1/npts);
  Fc=F*ctr;
end
