% Midterm exam II
% Manjul Regmi
% Q no. 1

clc; clear; close all;

a = 1;
b = 1;

p = 100;         % Penalty for tracking error at final time
q = 5;          % Penalty for tracking error at during
R = 2;
c = 1;

T = 20;
t = 0:0.1:T;
x0 = 3;

% Tracking signal
r = t.^2;
r_T = r(end);
% Boundary condition
s_T = c'*p*c;
v_T = c'*p*r_T;

[t1, s] = ode45(@(t,s)ricatti(t,s,a,b,c,q,R),t,s_T);

[t2, v] = ode45(@(t,v)aux(t,v,a,b,c,q,R,s,r,t1),t,v_T);

s = flip(s);
v = flip(v);

% Compute Kalman gains
K = zeros(1,length(t1));
for i = 1:length(t1)
    K(i) = inv(R)*b'*s(i);
end

[t3, x] = ode45(@(t,x)plant(t,x,a,b,R,K,v,t1),t,x0);

% Compute optimal control
u = zeros(1,length(t1));
for i = 1:length(t1)
    u(i) = -K(i)*x(i) + inv(R)*b'*v(i);
end

% Plot
figure
plot(t,r,'--')
legend('Reference signal')
figure
plot(t,s, 'linewidth', 2)
legend('Performance kernel')
figure
plot(t,K, 'linewidth', 2)
legend('Kalman gain')
figure
plot(t,v, 'linewidth', 2)
legend('Auxiliary signal')

figure
plot(t,x, 'linewidth', 2)
legend('State trajectory')

figure
plot(t,u, 'linewidth', 2)
legend('Control sequence')

% Functions
function [s_dot] = ricatti(t,s,a,b,c,q,R)
    s_dot = a'*s + s*a - s*b*inv(R)*b'*s +c'*q*c;
end

function [v_dot] = aux(t,v,a,b,c,q,R,s1,ref,t1)
    %Interpolate
    s = interp1(t1,s1,t);           % Performance kernel
    r = interp1(t1,ref,t);          % Reference trajectory

    K = inv(R)*b'*s;
    v_dot = (a-b*K)'*v + c'*q*r;    %r = reference tracking signal
end

function [x_dot] = plant(t,x,a,b,R,Kref,vref,t1)
    K = interp1(t1,Kref,t);
    v = interp1(t1,vref,t);

    u = -K*x + inv(R)*b'*v;
    x_dot = a*x + b*u;
end