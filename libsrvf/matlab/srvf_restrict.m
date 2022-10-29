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


% Returns a new SRVF representing the restriction of Q,T to the interval 
% [a,b].  [a,b] should be a subset of [T(1),T(end)], the domain of the SRVF.
%
% Inputs:
%  Q,T : the sample points and parameter values of the SRVF
%  a,b : the endpoints of the new domain interval
%
% Outputs:
%  Qr,Tr : the restricted SRVF
% --------------------------------------------------------------------------
function [Qr Tr] = srvf_restrict(Q,T,a,b)
  Tr = unique([T a b]);
  Qr = srvf_refine(Q,T,Tr);
  aidx = find(Tr==a);
  bidx = find(Tr==b);
  Qr = Qr(:,aidx:bidx-1);
  Tr = Tr(aidx:bidx);
end


%!test
%! T=[0 1/4 1/2 3/4 1];
%! Q=[0 1 0 1];
%! a=1/8;
%! b=5/8;
%! Trx=[1/8 1/4 1/2 5/8];
%! Qrx=[0 1 0];
%! [Qr Tr]=srvf_restrict(Q,T,a,b);
%! assert(Qr,Qrx);
%! assert(Tr,Trx);
