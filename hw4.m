clear 
clc

% Problem 12.1.4

% a.)

A = [10 -12 -6; 5 -5 -4; -1 0 3];
x = [1; 1; 1];
k = 10;
[lam1,u] = rqi(A,x,k);
x = [10; 1; 1];
[lam2,u] = rqi(A,x,k);
x = [1; 1; 33]; %krizia guessed this for me BOOM
[lam3,u] = rqi(A,x,k);
lam_a = [lam1; lam2; lam3];
 
% b.)

A = [-14 20 10; -19 27 12; 23 -32 -13];
x = [1; 1; 1];  
[lam1,u] = rqi(A,x,k);
x = [1; 1; 33];  
[lam2,u] = rqi(A,x,k);
x = [3; 0; 5];  
[lam3,u] = rqi(A,x,k);
lam_b = [lam1; lam2; lam3];
 
% c.)

A = [8 -8 -4; 12 -15 -7; -18 26 12];
x = [1; 0; 0];
[lam1,u] = rqi(A,x,k);
x = [1; 1; 1];
[lam2,u] = rqi(A,x,k);
x = [1; 10; 1];
[lam3,u] = rqi(A,x,k);
lam_c = [lam1; lam2; lam3];

% d.)

A = [12 -4 -2; 19 -19 -10; -35 52 27];
x = [1; 0; 0];
[lam1,u] = rqi(A,x,k);
x = [3; 1; 3];
[lam2,u] = rqi(A,x,k);
x = [1; 5; 2];
[lam3,u] = rqi(A,x,k);
lam_d = [lam1; lam2; lam3];

% Output for problem 12.1.4

lam_a % part a eigenvalues
lam_b % part b eigenvalues
lam_c % part c eigenvalues
lam_d % part d eigenvalues

% Problem 12.2.2

% a.)

a = [3 1 -2; 4 1 1; -3 0 3];
lam_a = shiftedqr0(a);

% b.)

a = [1 5 4; 2 -4 -3; 0 -2 4];
lam_b = shiftedqr0(a);

% c.)

a = [1 1 -2; 4 2 -3; 0 -2 2];
lam_c = shiftedqr(a);

% d.)

a = [5 -1 3; 0 6 1; 3 3 -3];
lam_d = shiftedqr0(a);

% Output for problem 12.2.2

lam_a % part a eigenvalues
lam_b % part b eigenvalues
lam_c % part c eigenvalues
lam_d % part d eigenvalues

% Problem 12.2.4

% a.)

a = [-1 1 3; 3 3 -2; -5 2 7];
[a1,v] = hessen(a);
lam_a = shiftedqr(a);

% b.)

a = [7 -33 -15; 2 26 7; -4 -50 -13];
[a2,v] = hessen(a);
lam_b = shiftedqr0(a);

% c.)

a = [8 0 5; -5 3 -5; 10 0 13];
[a3,v] = hessen(a);
lam_c = shiftedqr0(a);

% d.)

a = [-3 -1 1; 5 3 -1; -2 -2 0];
[a4,v] = hessen(a);
lam_d = shiftedqr(a);

% Output for problem 12.2.4

a1 % part a Hessenberg form matrix 
lam_a % part a eigenvalues
a2 % part b Hessenberg form matrix 
lam_b % part b eigenvalues
a3 % part c Hessenberg form matrix 
lam_c % part c eigenvalues
a4 % part d Hessenberg form matrix 
lam_d % part d eigenvalues

% Problem 12.4.4

% a.)

A = [1 2 2 1; 2 3 1 1; 4 5 6 3];
[U S V] = svd(A);
S1 = zeros(3,4);
S1(1) = S(1);
B1 = U*S1*V';
z = B1(3,:);
z = z';
x1 = B1(1:2,:);
x1 = x1';
x = ones(4,3);
x(:,1:2) = x1;
c1 = inv(x'*x)*(x'*z);

% b.)

A = [1 2 2 1; 2 3 1 1; 4 5 6 3];
[U S V] = svd(A);
S1 = zeros(3,4);
S1(1) = S(1);
B2 = U*S1*V';
z = B2(3,:);
z = z';
x1 = B2(1:2,:);
x1 = x1';
x = ones(4,3);
x(:,1:2) = x1;
c2 = inv(x'*x)*(x'*z);

% Output for problem 12.2.4

B1 % part a projected vectors
c1 % part a best least squares approximating plane
B2 % part a projected vectors
c2 % part a best least squares approximating plane

% Problem 12.4.6

% a.)

A1 = [3 0; 4 0];
[U1,S1,Vtrans1] = svd(A1); 

% b.)

A2 = [6 -2; 8 3/2];
[U2,S2,Vtrans2] = svd(A2); 

% c.)

A = [0 1; 0 0];

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

% d.)

A4 = [-4 -12; 12 11];
[U4,S4,Vtrans4] = svd(A4); 

% e.)

A5 = [0 -2; -1 0];
[U5,S5,Vtrans5] = svd(A5); 

% Output for problem 12.4.6

% a.)
U1 % U for part a
S1 % S for part a
Vtrans1 % V' for part a
U1*S1*Vtrans1 % proof that A = U*S*V'

% b.)
U2 % U for part b
S2 % S for part b
Vtrans2 % V' for part b
U2*S2*Vtrans2 % proof that A = U*S*V'

% c.)
U % U for part c
S % S for part c
Vtrans % V' for part c
U*S*Vtrans % proof that A = U*S*V'

% d.)
U4 % U for part d
S4 % S for part d
Vtrans4 % V' for part d
U4*S4*Vtrans4 % proof that A = U*S*V'

% e.)
U5 % U for part e
S5 % S for part e
Vtrans5 % V' for part e
U5*S5*Vtrans5 % proof that A = U*S*V'