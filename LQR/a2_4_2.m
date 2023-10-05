clc; clear; close all;

r = 0.1;
A = [1 1;0 1];
B = [0;1];
I = [1 0;0 1];
h = [1 -1 0 -10; 0 1 0 10; 1 -1 1 -10; 0 1 1 11];
[V,D] = eigs(h);

w11 = V(1:2,1:2);
w12 = V(1:2,3:4);
w21 = V(3:4,1:2);
w22 = V(3:4,3:4);

sn = [1 0;0 1];

t = -inv(w22-sn*w12) * (w21-sn*w11);

s_inf = w21 * inv(w11);

inv_h = inv(h);

[vv,dd] = eigs(inv_h);
M = dd(3:4,3:4);
X = vv(1:2,3:4);
lambda = vv(3:4,3:4);

K_inf = (1/r)*B'*lambda*M*inv(X);

% Ackermann formula
U2 = [B A*B];
deltaA = A*A - (trace(M))*A + M(1,1)*M(2,2)*I;
K = [0 1]*inv(U2)*deltaA;