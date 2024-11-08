clear all
close all
clc

% Define symbolic variables
syms x0(t) x1(t) x2(t) x3(t) y3(t) th3(t) ...
    k1 k2 k3 kt3 k4 k5 ...
    c1 c2 c3 ct3 c4 c5 ...
    m1 m2 m3 J3g ...
    L Lg th30 ...
    omega2 lambda

% Differentiate (state) variables
x0dot = diff(x0,t);
x1dot = diff(x1,t);
x2dot = diff(x2,t);
x3dot = diff(x3,t);
y3dot = diff(y3,t);
th3dot = diff(th3,t);

% Put state variables in state vector q
q = [x1 x2 x3 y3 th3];

% Put differentiated state variables in  vector qdot
qdot = [x1dot x2dot x3dot y3dot th3dot];


% Potential energy
V = 1/2*(k1*(x1(t)-x0(t))^2 + k2*(x2(t)-x1(t))^2 + k3*(x3(t)-x2(t))^2 +kt3*(th3(t)-th30)^2 ...
    + k4*(x3(t)+L*cos(th3(t))-L*cos(th30))^2 +k5*(y3(t)+L*sin(th3(t))-L*sin(th30))^2);

% Kinetic energy
T = 1/2*(m1*x1dot^2 + m2*x2dot^2 + m3*(x3dot^2 + y3dot^2 + 2*Lg*y3dot*th3dot*cos(th3) - ...
    2*Lg*x3dot*th3dot*sin(th3) + Lg^2*th3dot^2) + J3g*th3dot^2);

% Dissipation function
D = 1/2*(c1*(x1dot-x0dot)^2 + c2*(x2dot-x1dot)^2 + c3*(x3dot-x2dot)^2 + ct3*th3dot^2 ...
    + c4*(x3dot-L*th3dot*sin(th3))^2 + c5*(y3dot+L*th3dot*cos(th3))^2);

% Lagrange

% Variable: x1
Tdot_x1dot = diff(T,x1dot);
Tdot_x1dot_t = diff(Tdot_x1dot,t);
Tdot_x1 = diff(T,x1);
Vdot_x1 = diff(V,x1);
Ddot_x1dot = diff(D,x1dot);

lagr_x1 = Tdot_x1dot_t - Tdot_x1 + Vdot_x1 + Ddot_x1dot;
% pretty(lagr_x1)

% Variable: x2
Tdot_x2dot = diff(T,x2dot);
Tdot_x2dot_t = diff(Tdot_x2dot,t);
Tdot_x2 = diff(T,x2);
Vdot_x2 = diff(V,x2);
Ddot_x2dot = diff(D,x2dot);

lagr_x2 = Tdot_x2dot_t - Tdot_x2 + Vdot_x2 + Ddot_x2dot;
% pretty(lagr_x2)

% Variable: x3
Tdot_x3dot = diff(T,x3dot);
Tdot_x3dot_t = diff(Tdot_x3dot,t);
Tdot_x3 = diff(T,x3);
Vdot_x3 = diff(V,x3);
Ddot_x3dot = diff(D,x3dot);

lagr_x3 = Tdot_x3dot_t - Tdot_x3 + Vdot_x3 + Ddot_x3dot;
% pretty(lagr_x3)

% Variable: y3
Tdot_y3dot = diff(T,y3dot);
Tdot_y3dot_t = diff(Tdot_y3dot,t);
Tdot_y3 = diff(T,y3);
Vdot_y3 = diff(V,y3);
Ddot_y3dot = diff(D,y3dot);

lagr_y3 = Tdot_y3dot_t - Tdot_y3 + Vdot_y3 + Ddot_y3dot;
% pretty(lagr_y3)

% Variable: th3
Tdot_th3dot = diff(T,th3dot);
Tdot_th3dot_t = diff(Tdot_th3dot,t);
Tdot_th3 = diff(T,th3);
Vdot_th3 = diff(V,th3);
Ddot_th3dot = diff(D,th3dot);

lagr_th3 = Tdot_th3dot_t - Tdot_th3 + Vdot_th3 + Ddot_th3dot;
% pretty(lagr_th3)


% Making the linearised matrices
Kvec = jacobian(V,q);
K = jacobian(Kvec,q);
% K2 = hessian(V,q);

Cvec = jacobian(D,qdot);
C = jacobian(Cvec,qdot);
% C2 = hessian(D,qdot);

Mvec = jacobian(T,qdot);
M = jacobian(Mvec,qdot);
% M2 = hessian(T,qdot);

% Substituting equilibrium values in matrices

% Value variables 
m1_v = 0.45;
m2_v = 1.15;
m3_v = 1.9;
J3g_v = 0.0149;
k1_v = 155.8;
k2_v = 23.6;
k3_v = 444.6;
k4_v = 415.4;
k5_v = 50.25;
kt3_v =2; 
c1_v = 30;
c2_v = 202.84;
c3_v = 500;
c4_v = 164.59;
c5_v = 49.99;
ct3_v = 4.9;
L_v = 0.298;
Lg_v = 0.6*L_v;
th30_v = pi/2;

% Value state variables around equilibrium
x1_v = 0;
x2_v = 0;
x3_v = 0;
y3_v = 0;
th3_v = pi/2;

Ks = subs(K,{k1,k2,k3,k4,k5,kt3,L,x3,y3,th3,th30},{k1_v,k2_v,k3_v,k4_v,k5_v,kt3_v,L_v,x3_v,y3_v,th3_v,th30_v});

Cs = subs(C,{c1,c2,c3,c4,c5,ct3,L,th3,th30},{c1_v,c2_v,c3_v,c4_v,c5_v,ct3_v,L_v,th3_v,th30_v});

Ms = subs(M,{m1,m2,m3,Lg,L,J3g,th3},{m1_v,m2_v,m3_v,Lg_v,L_v,J3g_v,th3_v});


% Eigenmodes and -frequencies without damping
poly = det(Ks-omega2.*Ms);
omega2_v = double(solve(poly,omega2)); % With vpa 'omega2' will stay symbolic. Use 'double' to make it an integer

% poly2 = det(lambda^2*Ms+Ks);
% lambda_v = double(solve(poly2,lambda));
% lambda_v = lambda_v.^2;

mat1 = double(Ks-omega2_v(1)*Ms); % Symbolic matrix to 'double' matrix
mat2 = double(Ks-omega2_v(2)*Ms);
mat3 = double(Ks-omega2_v(3)*Ms);
mat4 = double(Ks-omega2_v(4)*Ms);
mat5 = double(Ks-omega2_v(5)*Ms);

mat1_rref = rref(mat1);
X1 = [-mat1_rref(1:4,end);1]/-mat1_rref(1,end);

mat2_rref = rref(mat2);
X2 = [-mat2_rref(1:4,end);1]/-mat2_rref(1,end);

mat3_rref = rref(mat3);
X3 = [-mat3_rref(1:4,end);1]/-mat3_rref(1,end);

mat4_rref = rref(mat4);
X4 = [-mat4_rref(1:4,end);1]/-mat4_rref(1,end);

mat5_rref = rref(mat5);
X5 = [-mat5_rref(1:4,end);1]/-mat5_rref(1,end);

X = [X1 X2 X3 X4 X5];
freq = diag(sqrt(omega2_v));

% Eigenmodes and -frequencies with damping
A = [zeros(5,5) Ms;Ms Cs];
B = [-Ms zeros(5,5); zeros(5,5) Ks];
poly_d = det(lambda*A + B);
lambda_v = double(solve(poly_d,lambda));

mat1d = double(Ks + lambda_v(1)*Cs + (lambda_v(1)^2)*Ms);
mat2d = double(Ks + lambda_v(2)*Cs + (lambda_v(2)^2)*Ms);
mat3d = double(Ks + lambda_v(3)*Cs + (lambda_v(3)^2)*Ms);
mat4d = double(Ks + lambda_v(4)*Cs + (lambda_v(4)^2)*Ms);
mat5d = double(Ks + lambda_v(5)*Cs + (lambda_v(5)^2)*Ms);
mat6d = double(Ks + lambda_v(6)*Cs + (lambda_v(6)^2)*Ms);
mat7d = double(Ks + lambda_v(7)*Cs + (lambda_v(7)^2)*Ms);
mat8d = double(Ks + lambda_v(8)*Cs + (lambda_v(8)^2)*Ms);
mat9d = double(Ks + lambda_v(9)*Cs + (lambda_v(9)^2)*Ms);
mat10d = double(Ks + lambda_v(10)*Cs + (lambda_v(10)^2)*Ms);

mat1d_rref = rref(mat1d);
X1d = [-mat1d_rref(1:4,end);1]/-mat1d_rref(1,end);

mat2d_rref = rref(mat2d);
X2d = [-mat2d_rref(1:4,end);1]/-mat2d_rref(1,end);

mat3d_rref = rref(mat3d);
X3d = [-mat3d_rref(1:4,end);1]/-mat3d_rref(1,end);

mat4d_rref = rref(mat3d);
X4d = [-mat4d_rref(1:4,end);1]/-mat4d_rref(1,end);

mat5d_rref = rref(mat5d);
X5d = [-mat5d_rref(1:4,end);1]/-mat5d_rref(1,end);

mat6d_rref = rref(mat6d);
X6d = [-mat6d_rref(1:4,end);1]/-mat6d_rref(1,end);

mat7d_rref = rref(mat7d);
X7d = [-mat7d_rref(1:4,end);1]/-mat7d_rref(1,end);

mat8d_rref = rref(mat8d);
X8d = [-mat8d_rref(1:4,end);1]/-mat8d_rref(1,end);

mat9d_rref = rref(mat9d);
X9d = [-mat9d_rref(1:4,end);1]/-mat9d_rref(1,end);

mat10d_rref = rref(mat10d);
X10d = [-mat10d_rref(1:4,end);1]/-mat10d_rref(1,end);

Xd = [X1d X2d X3d X4d X5d X6d X7d X8d X9d X10d];

freq_d = sqrt(-lambda_v);

% Modal masses for each eigenfrequency and -mode
mmass1 = double(X1d.'*Ms*X1d);
mmass2 = double(X2d.'*Ms*X2d);
mmass3 = double(X3d.'*Ms*X3d);
mmass4 = double(X4d.'*Ms*X4d);
mmass5 = double(X5d.'*Ms*X5d);
mmass6 = double(X6d.'*Ms*X6d);
mmass7 = double(X7d.'*Ms*X7d);
mmass8 = double(X8d.'*Ms*X8d);
mmass9 = double(X9d.'*Ms*X9d);
mmass10 = double(X10d.'*Ms*X10d);

% Modal damping factor for each eigenfrequency and -mode
Beta1 = double(X1d.'*Cs*X1d);
Beta2 = double(X2d.'*Cs*X2d);
Beta3 = double(X3d.'*Cs*X3d);
Beta4 = double(X4d.'*Cs*X4d);
Beta5 = double(X5d.'*Cs*X5d);
Beta6 = double(X6d.'*Cs*X6d);
Beta7 = double(X7d.'*Cs*X7d);
Beta8 = double(X8d.'*Cs*X8d);
Beta9 = double(X9d.'*Cs*X9d);
Beta10 = double(X10d.'*Cs*X10d);

%  Modal damping ratios for eacht eigenfrequency and - mode
Epsilon1 = Beta1/(2*freq_d(1)*mmass1);
Epsilon2 = Beta2/(2*freq_d(2)*mmass2);
Epsilon3 = Beta3/(2*freq_d(3)*mmass3);
Epsilon4 = Beta4/(2*freq_d(4)*mmass4);
Epsilon5 = Beta5/(2*freq_d(5)*mmass5);
Epsilon6 = Beta6/(2*freq_d(6)*mmass6);
Epsilon7 = Beta7/(2*freq_d(7)*mmass7);
Epsilon8 = Beta8/(2*freq_d(8)*mmass8);
Epsilon9 = Beta9/(2*freq_d(9)*mmass9);
Epsilon10 = Beta10/(2*freq_d(10)*mmass10);



% Generate motion 
damping_type = 2; % 1 for Undamped, 2 for damped.
chosen_mode = 10; % Selected eigen mode
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
x4 = x3 + disturb .* (L_v * cos(eigenvector(5))) + initial(5);
y4 = y3 + disturb .* (L_v * sin(eigenvector(5))) + L_v * 100;

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
filename_prefix = "EDassignment41"; % Adjust prefix
filename = join([filename_prefix, damping_type_str, mode_number_str]);
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




