% libsrvf 
% =======
%
% A shape analysis library using the square root velocity framework.
% 
% Copyright (C) 2014   FSU Statistical Shape Analysis and Modeling Group
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


% Returns a reparametrization which yields a constant-speed parametrized 
% function when applied to the input PLF.
%
% Inputs:
%  F,T - the sample points and parameter values of the PLF.
%        F is a matrix containing one column per sample point.
% 
% Outputs:
%  G,T - the sample points and parameter values of the reparametrization
% --------------------------------------------------------------------------
function [G GT] = plf_constant_speed_reparam(F, T)
  dF = diff(F, 1, 2);
  GT = [0 cumsum(sqrt(sum(dF .* dF, 1)))];
  GT = GT ./ GT(end);
  G = T;
end
