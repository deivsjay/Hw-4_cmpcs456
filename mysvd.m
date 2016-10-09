function [U,S,Vtrans] = mysvd(A) 

[U,D] = eig(A*A');
temp = U(:,2);
U(:,2) = U(:,1);
U(:,1) = temp;

[V,D] = eig(A'*A);
temp = V(:,2);
V(:,2) = V(:,1);
V(:,1) = temp;
Vtrans = V';

e = eig(A*A');
temp = e(2,:);
e(2,:) = e(1,:);
e(1,:) = temp;

S = zeros(2,2);
S(1,1) = sqrt(e(1));
S(2,2) = sqrt(e(2));

