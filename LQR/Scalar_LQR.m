% Close loop feedback control with free final state
% LQR control
% Scalar case

N = 5;     % Final time

S = zeros(1,N+1);
S(N+1) = 10;   % Penalty for final state
K = zeros(1,N);


a = 2;
b = 1;
q = 1;
r = 1;

for i=N:-1:1
    S(i) = (a^2*r*S(i+1))/(r + b^2*S(i+1)) + q;
    K(i) = (a*b*S(i+1))/(r + b^2*S(i+1));
end

x = zeros(1,N+1);
u = zeros(1,N);
x(1) = 3;     % Initial state
for i=1:N
    u(i) = -K(i)*x(i);
    x(i+1) = a*x(i) + b*u(i);
end
