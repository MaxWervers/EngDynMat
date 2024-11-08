clear;
clc;
syms m1 m2 m3 L Lg J3g x0(t) x1(t) x2(t) x3(t) y3(t) theta3(t) theta3_hat t k1 k2 k3 k4 k5 kt3 c1 c2 c3 ct3 c4 c5


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
          (1/2)*m3*(x3_dot - Lg*theta3_dot*sin(theta3))^2 + ...
          (1/2)*m3*(y3_dot + Lg*theta3_dot*cos(theta3))^2 + ...
          (1/2)*(J3g)*theta3_dot^2;
         
% Define the total potential energy expression V_total
V_total = (1/2)*k1*(x0 - x1)^2 + ...
          (1/2)*k2*(x1 - x2)^2 + ...
          (1/2)*k3*(x2 - x3)^2 + ...
          (1/2)*kt3*(theta3 - theta3_hat)^2 + ...
          (1/2)*k4*(x3 + L*cos(theta3) - (L*cos(theta3_hat)))^2 + ...
          (1/2)*k5*(y3 + L*sin(theta3) - (L*sin(theta3_hat)))^2;
          
% Define the total dissipation expression D_total
D_total = (1/2)*c1*(x0_dot - x1_dot)^2 + ...
          (1/2)*c2*(x1_dot - x2_dot)^2 + ...
          (1/2)*c3*(x2_dot - x3_dot)^2 + ...
          (1/2)*ct3*(theta3_dot)^2 + ...
          (1/2)*c4*(x3_dot - L*theta3_dot*sin(theta3))^2 + ...
          (1/2)*c5*(y3_dot + L*theta3_dot*cos(theta3))^2;

% Langrangian eqations of motion:
L_x1 = diff(T_total, x1_dot, t) - diff(T_total, x1) + diff(V_total, x1) + diff(D_total, x1_dot);
L_x2 = diff(T_total, x2_dot, t) - diff(T_total, x2) + diff(V_total, x2) + diff(D_total, x2_dot);
L_x3 = diff(T_total, x3_dot, t) - diff(T_total, x3) + diff(V_total, x3) + diff(D_total, x3_dot);
L_y3 = diff(T_total, y3_dot, t) - diff(T_total, y3) + diff(V_total, y3) + diff(D_total, y3_dot);
L_theta3 = diff(T_total, theta3_dot, t) - diff(T_total, theta3) + diff(V_total, theta3) + diff(D_total, theta3_dot);

% Display the Lagrangian equations of motion
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

% Define the generalized coordinates and their derivatives
q_vec = [x1 ,x2, x3, y3, theta3];
q_dot_vec = [x1_dot, x2_dot, x3_dot, y3_dot, theta3_dot];

% Define the mass matrix M, the Damping matrix C, and the stiffness matrix K
M = hessian(T_total, q_dot_vec)
C = hessian(D_total, q_dot_vec)
K = hessian(V_total, q_vec)

% Define the mass matrix M, the Damping matrix C, and the stiffness matrix K at q_eq = [0, 0, 0, 0, pi/2]
M_at_theta3 = subs(M, theta3, pi/2)

C_at_theta3 = subs(C, theta3, pi/2)

K_at_theta3 = subs(subs(subs(K, theta3, pi/2), y3, 0), theta3_hat, pi/2)

% Define the generalized forces Q
Q = mtimes(C_at_theta3, transpose([x0_dot ,0 ,0 ,0 ,0])) + mtimes(K_at_theta3, transpose([x0 ,0 ,0 ,0 ,0])) 

% Call the function to get the model parameters
params = getModelParameters();

% Access the necessary values
m1_num = params.m1;
m2_num = params.m2;
m3_num = params.m3;
Lg_num = params.Lg;
J3g_num = params.J3g;

k1_num = params.k1;
k2_num = params.k2;
k3_num = params.k3;
kt3_num = params.kt3;
k4_num = params.k4;
k5_num = params.k5;
L_num = params.L;

c1_num = params.c1;
c2_num = params.c2;
c3_num = params.c3;
ct3_num = params.ct3;
c4_num = params.c4;
c5_num = params.c5;





% Substitute all variables into K, C and M matrices
M_substituted = subs(M_at_theta3, {m1, m2, m3, Lg, J3g}, {m1_num, m2_num, m3_num, Lg_num, J3g_num})

K_substituted = subs(K_at_theta3, {k1, k2, k3, kt3, k4, k5, L}, {k1_num, k2_num, k3_num, kt3_num, k4_num, k5_num, L_num})

C_substituted = subs(C_at_theta3, {c1, c2, c3, ct3, c4, c5, L}, {c1_num, c2_num, c3_num, ct3_num, c4_num, c5_num, L_num})

% Conversion to numeric
M_num = double(M_substituted)
K_num = double(K_substituted)
C_num = double(C_substituted)

%Find eigenfrequencies & eigenmodes of undamped system
[eigenmodes_undamped, omega_squared] = eig(K_num, M_num)

eigen_freq_undamped = sqrt(omega_squared)

%Find eigenfrequencies & eigenmodes of damped system
[eigenmodes_damped, eigen_freq_damped] = polyeig(K_num,C_num,M_num)

eigen_freq_damped = sqrt(eigen_freq_damped)

%find modal mass matrix
mu = diag(transpose(eigenmodes_damped) * M_num * eigenmodes_damped)

%find modal damping matrix
beta = diag(transpose(eigenmodes_damped) * C_num * eigenmodes_damped)

%find modal stiffness matrix
gamma = diag(transpose(eigenmodes_damped) * K_num * eigenmodes_damped)

%find damping ratios
epsilon = beta ./ (2.*eigen_freq_damped.*mu)

% initialize for animation
X = eigenmodes_undamped;
Xd = eigenmodes_damped;

omega2_v = diag(eigen_freq_undamped)
lambda_v = eigen_freq_damped


% Generate motion 
damping_type = 1; % 1 for Undamped, 2 for damped.
chosen_mode = 1; % Selected eigen mode
time_values = 0:0.01:1; % Time values

% Introduce disturbance
initial = [0, 10, 25, 0, 0]; 

% Damping yes/no
if damping_type == 1
    damping_type_str = "Undamped";
    eigenvector = real(X(:, chosen_mode));
    disturb = 0.5 * sin(abs(real(omega2_v(chosen_mode))) * time_values);
    title_text = "Mode " + chosen_mode + " without damping for \omega = " ...
        + double(omega2_v(chosen_mode)) + " rad/s";
else
    damping_type_str = "Damped";
    eigenvector = real(Xd(:, chosen_mode));
    disturb = 0.5 * sin(abs(real(lambda_v(chosen_mode))) * time_values);
    title_text = "Mode " + chosen_mode + " with damping for \omega = " ...
        + double(lambda_v(chosen_mode)) + " rad/s";
end

% Integrate the disturbed system
x1 = eigenvector(1) * disturb + initial(1);
x2 = eigenvector(2) * disturb + initial(2);
x3 = eigenvector(3) * disturb + initial(3);
y3 = eigenvector(4) * disturb + initial(4);
x4 = x3 + disturb .* (L_num * cos(eigenvector(5))) + initial(5);
y4 = y3 + disturb .* (L_num * sin(eigenvector(5))) + L_num * 100;

% Draw initial figure
mass1 = plot([x1(1) x2(1)], [0 0], 'o', 'MarkerSize', 5);
hold on
mass2 = plot([x3(1) x4(1)], [y3(1) y4(1)], 'o', 'MarkerSize', 5); 
hold on
xlim([-20, 50]);
ylim([-20, 50]);
xlabel('X(mm)') 
ylabel('Y(mm)')
title(title_text)

% Animation loop
mode_number_str = ".Mode" + chosen_mode; % Adjust label
filename_prefix = "EDassignment10"; % Adjust prefix
filename = join(["Simulations/", filename_prefix, damping_type_str, mode_number_str]);
animation_video = VideoWriter(filename, 'MPEG-4');
animation_video.FrameRate = 10;
open(animation_video)

for frame_index = 1:length(disturb)
    % Update mass positions
    set(mass1, 'XData', [x1(frame_index) x2(frame_index)], 'YData', [0 0]);
    hold on
    set(mass2, 'XData', [x3(frame_index) x4(frame_index)], 'YData', [y3(frame_index) y4(frame_index)]);
    
    % Draw lower
    lower = line([x1(frame_index), x2(frame_index), x3(frame_index)], [0, 0, y3(frame_index)], ...
                           'Color', [0.8 0.2 0.6], 'Linewidth', 1); 
                           
    % Draw upper
    upper = line([x3(frame_index), x4(frame_index)], [y3(frame_index), y4(frame_index)], ...
                          'Color', [0.2 0.6 0.8], 'Linewidth', 2.5);
    
    % Capture the frame
    frame = getframe(gcf);
    image_data = frame2im(frame);
    writeVideo(animation_video, image_data);
    
    % Remove drawn elements to prepare for the next frame
    delete(lower)
    delete(upper)
end

close(animation_video)