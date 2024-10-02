function params = getModelParameters()
    % Function to return the arm model parameters in a structured format

    % Mass parameters
    params.m1 = 0.5;       % kg, Mass of palm
    params.m2 = 1.35;      % kg, Mass of lower arm
    params.m3 = 2.1;       % kg, Mass of upper arm

    % Mass moment of inertia
    params.J3g = 0.016;    % kg*m^2, Mass moment of inertia of upper arm

    % Spring constants
    params.k1 = 195.8 * 1e3;  % N/m, Spring constant of device-palm interactions (converted from kN/m)
    params.k2 = 23.61 * 1e3;  % N/m, Spring constant of palm-lower arm interactions (converted from kN/m)
    params.k3 = 484.6 * 1e3;  % N/m, Spring constant of lower arm-elbow interactions (converted from kN/m)
    params.k4 = 515.4 * 1e3;  % N/m, Spring constant of upper arm-shoulder interactions in x-direction (converted from kN/m)
    params.k5 = 40.25 * 1e3;  % N/m, Spring constant of upper arm-shoulder interactions in y-direction (converted from kN/m)

    % Rotational spring constant
    params.kt3 = 2.100;    % N*m/rad, Rotational stiffness of elbow

    % Viscous damping constants
    params.c1 = 32.45;     % N*s/m, Viscous damping of device-palm interactions
    params.c2 = 160.8;     % N*s/m, Viscous damping of palm-lower arm interactions
    params.c3 = 593.6;     % N*s/m, Viscous damping of lower arm-elbow interactions
    params.c4 = 134.6;     % N*s/m, Viscous damping of upper arm-shoulder interactions in x-direction
    params.c5 = 50.01;     % N*s/m, Viscous damping of upper arm-shoulder interactions in y-direction

    % Rotational viscous damping constant
    params.ct3 = 4.726;    % N*m*s/rad, Rotational viscous damping of elbow

    % Lengths
    params.Lg = 0.55 * 0.328; % m, Distance between elbow and mass center of upper arm
    params.L = 0.328;         % m, Length of upper arm

    % Operational angel
    params.theta3_hat = pi/2;      % rad, To be defined
end
