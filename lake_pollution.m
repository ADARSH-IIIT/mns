clc; clear; close all;

%% Parameters
flowin = 100;        % Inflow rate (m^3/s)
flowout = 100;       % Outflow rate (m^3/s) (steady-state condition)
Cin = 5;            % Pollutant concentration of inflow (mg/L)
V = 5000;           % Volume of lake (m^3)
t0 = 2;            % Initial time (not zero)
C0 = 10;           % Initial concentration of pollution in lake (mg/L)
tf = 150;           % Final simulation time
h = 1;             % Step size for Euler's method

%% Euler's Method
t_euler = t0:h:tf;  % Time vector
C_euler = zeros(size(t_euler));
C_euler(1) = C0; % Initial condition

% Euler Iteration
for i = 2:length(t_euler)
    dCdt = (flowin * Cin - flowout * C_euler(i-1)) / V;
    C_euler(i) = C_euler(i-1) + h * dCdt;
end

%% Using ode45
% Define the function for ode45
lake_pollution_fun = @(t, C) (flowin * Cin - flowout * C) / V;
tspan = [t0 tf];

% Solve ODE using ode45
[t_ode, C_ode] = ode45(lake_pollution_fun, tspan, C0);

%% Plot Results
figure;
hold on;
plot(t_euler, C_euler, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(t_ode, C_ode, 'b-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Pollutant Concentration (mg/L)');
title('Lake Pollution Model: Euler vs. ode45');
legend('Euler Method', 'ode45');
grid on;
hold off;
