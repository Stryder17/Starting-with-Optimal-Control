% Homework 5
% Manjul Regmi(2254052)

clc; clear; close all;

% System Dynamics from 3.4-1
A = [0 1; 0 0];
B = [0; 1];
C = [1 0];
D= 0;

% Problem 4.2-2
nu = 0;
q = 10;
Q = [1 nu; nu q];
R = 1;
T = 5;
x0 = [0; 0];
%x_T = 10;
r_T = 10;

% Boundary conditions
S_T = [0 0;0 0];
V_T = C';
P_T = 0;

K0 = [0 0];      % Set Kalman gain to zero initially
% Convert boundary conditions to vector to use in ode45
boundary = [S_T(1,1); S_T(1,2); S_T(2,2); V_T; P_T; K0'];

tspan = 0:0.02:T;
% Integrate using ode45
[t1,var] = ode45(@(t,var)optimalControl(t,var,A,B,Q,R), tspan, boundary);
%V is a vector of dimension n(Column)
%P is a vector. For 1 control input this vector becomes scalar.

% Solution of differential equations
% The solution varies from final time to initial time. So flip the elements.
S(1,1,:) = flip(var(:,1));
S(1,2,:) = flip(var(:,2));
S(2,1,:) = flip(var(:,2));
S(2,2,:) = flip(var(:,3));

V(1,:) = flip(var(:,4));
V(2,:) = flip(var(:,5));

P(1,:) = flip(var(:,6));

K(:,1) = flip(var(:,7));
K(:,2) = flip(var(:,8));

%t = flip(t1);
t = t1;

% Obtain control sequence and state trajectory
u = zeros(1,length(t1));
x = zeros(2,length(t1));
x(:,1) = x0;                            % Initial condition

% for i=1:length(t1)-1
%     u(i) = -(K(i,:)-inv(R)*B'*V(:,i)*inv(P(i))*V(:,i)')*x(:,i) - inv(R)*B'*V(:,i)*inv(P(i))*r_T;
%     del_T = t(i+1)-t(i);
%     x(:,i+1) = x(:,1) + (A*x(:,i) + B*u(i))*del_T;          % Using finite difference method
% 
% end

for i=1:length(t1)-1
    u(i) = -(K(i,:)-inv(R)*B'*V(:,i)*inv(P(i))*V(:,i)')*x(:,i) - inv(R)*B'*V(:,i)*inv(P(i))*r_T;
    del_T = t(i+1)-t(i);
    x(:,i+1) = expm(A - B*(K(i,:)-inv(R)*B'*V(:,i)*inv(P(i))*V(:,i)'))*del_T*x(:,i) - (B* inv(R)*B'*V(:,i)*inv(P(i))*r_T)*del_T;

end


% Plot control sequence and state trajectory
figure
plot(t,u)

figure
plot(t,x(1,:))

figure
plot(t,x(2,:))

% sys=ss(A,B,C,D);
% [X,t]=lsim(sys,u,t,x0); %Simulation with the optimal control
% figure
% plot(t,X)
% %-----------------------------------------------------------------------%
function [vec] = optimalControl(t,var,A,B,Q,R)
    % Unpack the vector
    S = [var(1) var(2); var(2) var(3)];
    V = [var(4); var(5)];
    P = var(6);
    
    % Differential equations
    SDot=+(A'*S + S*A - S*B*inv(R)*B'*S + Q);
    K = inv(R)*B'*S;
    VDot = +(A-B*K)'*V;
    PDot = -(V'*B*inv(R)*B'*V);

    % Collect all values into a single vector
    vec = [SDot(1,1); SDot(1,2); SDot(2,2); VDot; PDot; K'];
end

function [x_dot] = system(t,x,A,B,K,R,V,P,r_T)
    x_dot = A*x - B*(K-inv(R)*B'*V*(1/P)*V')*x - B* inv(R)*B'*V*(1/P)*r_T;
end