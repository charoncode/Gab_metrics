function v=curveDeriv(c)

% Input: plane curve c as 2xn matrix
% Output: discrete derivative curve

[d,n]=size(c);

for i=1:d
    v(i,:)=gradient(c(i,:),1/n);
end