clear;
clc;
syms m1 m2 m3 L Lg J3g x0(t) x1(t) x2(t) x3(t) y3(t) theta3(t) theta3_hat t k1 k2 k3 k4 k5 kt3 c1 c2 c3 ct3 c4 c5

%%%%%% Im not certain how we should define x0_dot %%%%%%
% Define symbolic derivatives with respect to time t
x0_dot = diff(x0, t);            % x0_dot = dx0/dt
x1_dot = diff(x1, t);            % x1_dot = dx1/dt
x2_dot = diff(x2, t);            % x2_dot = dx2/dt
x3_dot = diff(x3, t);            % x3_dot = dx3/dt
y3_dot = diff(y3, t);            % y3_dot = dy3/dt
theta3_dot = diff(theta3, t);    % theta3_dot = dtheta3/dt


% Define the total kinetic energy expression T_total
T_total = (1/2)*m1*x1_dot^2 + ...
          (1/2)*m2*x2_dot^2 + ...
          (1/2)*m3*((Lg/L)*x3_dot + (L - Lg)/L * (x3_dot - L*theta3_dot*sin(theta3)))^2 + ...
          (1/2)*m3*((Lg/L)*y3_dot + (L - Lg)/L * (theta3_dot + L*theta3_dot*cos(theta3)))^2 + ...
          (1/2)*J3g*theta3_dot^2;

% Define the total potential energy expression V_total
V_total = (1/2)*k1*(x0 - x1)^2 + ...
          (1/2)*k2*(x1 - x2)^2 + ...
          (1/2)*k3*(x2 - x3)^2 + ...
          (1/2)*kt3*(theta3 - theta3_hat)^2 + ...
          (1/2)*k4*(x3 + L*cos(theta3))^2 + ...
          (1/2)*k5*(y3 + L*sin(theta3))^2;
          
% Define the total dissipation expression D_total
D_total = (1/2)*c1*(x0_dot - x1_dot)^2 + ...
          (1/2)*c2*(x1_dot - x2_dot)^2 + ...
          (1/2)*c3*(x2_dot - x3_dot)^2 + ...
          (1/2)*ct3*(theta3_dot)^2 + ...
          (1/2)*c4*(y3_dot + L*theta3_dot*cos(theta3))^2 + ...
          (1/2)*c5*(x3_dot + L*theta3_dot*sin(theta3))^2

% Langrangian eqations:
L_x1 = diff(T_total, x1_dot, t) - diff(T_total, x1) + diff(V_total, x1) + diff(D_total, x1);
L_x2 = diff(T_total, x2_dot, t) - diff(T_total, x2) + diff(V_total, x2) + diff(D_total, x2);
L_x3 = diff(T_total, x3_dot, t) - diff(T_total, x3) + diff(V_total, x3) + diff(D_total, x3);
L_y3 = diff(T_total, y3_dot, t) - diff(T_total, y3) + diff(V_total, y3) + diff(D_total, y3);
L_theta3 = diff(T_total, theta3_dot, t) - diff(T_total, theta3) + diff(V_total, theta3) + diff(D_total, theta3);

% Display the Lagrangian equations
disp('The Lagrangian equation for x1 is:');
disp(L_x1);

disp('The Lagrangian equation for x2 is:');
disp(L_x2);

disp('The Lagrangian equation for x3 is:');
disp(L_x3);

disp('The Lagrangian equation for y3 is:');
disp(L_y3);

disp('The Lagrangian equation for theta3 is:');
disp(L_theta3);


