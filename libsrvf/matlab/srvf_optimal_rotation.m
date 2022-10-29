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


% Given two SRVFs Q1 and Q2, returns the rotation R which best aligns Q2 to Q1.
%
% Inputs:
%  Q1,T1 : the sample points and parameter values of the first SRVF
%  Q2,T2 : the sample points and parameter values of the second SRVF
%
% Outputs:
%  R : a D x D matrix (where D is the dimension of Q1 and Q2) representing 
%      the rotation.
% --------------------------------------------------------------------------
function R = srvf_optimal_rotation( Q1, T1, Q2, T2 )
  dim = size(Q1,1);

  % Rotations don't make sense for 1-D functions
  if (dim < 2)
    R = [1];
    return;
  end

  Tr = unique( [T1 T2] );
  Q1r = srvf_refine( Q1, T1, Tr );
  Q2r = srvf_refine( Q2, T2, Tr );

  for i=1:dim
    Q2r(i,:) = Q2r(i,:) .* diff( Tr );
  end

  A = Q1r * Q2r';
  [U,S,V] = svd(A);
  if det(A) > 0
    S = eye(dim);
  else
    S = eye(dim);
    S(:,end) = -S(:,end);
  end
  R = U*S*V';
end


%!demo
%! load ../tests/data/horse-1.csv
%! load ../tests/data/horse-2.csv
%! X1=horse_1;
%! X2=horse_2;
%! T1 = linspace(0,1,length(X1));
%! T2 = linspace(0,1,length(X2));
%! Q1=plf_to_srvf(X1,T1);
%! Q2=plf_to_srvf(X2,T2);
%! R=srvf_optimal_rotation(Q1,T1,Q2,T2);
%! X2r=R*X2;
%!
%! figure();
%! plot(X1(1,:),X1(2,:),'b',X2(1,:),X2(2,:),'r');
%! xl=min([X1(1,:) X2(1,:)]);
%! xu=max([X1(1,:) X2(1,:)]);
%! yl=min([X1(2,:) X2(2,:)]);
%! yu=max([X1(2,:) X2(2,:)]);
%! axis([xl xu yl yu],'square');
%! title('Curves before rotational alignment');
%!
%! figure();
%! plot(X1(1,:),X1(2,:),'b',X2r(1,:),X2r(2,:),'r');
%! xl=min([X1(1,:) X2r(1,:)]);
%! xu=max([X1(1,:) X2r(1,:)]);
%! yl=min([X1(2,:) X2r(2,:)]);
%! yu=max([X1(2,:) X2r(2,:)]);
%! axis([xl xu yl yu],'square');
%! title('Curves after rotational alignment');
