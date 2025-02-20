clc; clear; close all;

%% Parameters
r = 0.5;       % Growth rate
K = 100;       % Carrying capacity
t0 = 2;       % Initial time (not zero)
P0 = 10;      % Initial population
tf = 20;      % Final time
h = 0.5;      % Step size for Euler's method

%% Euler's Method
t_euler = t0:h:tf;  % Time vector
P_euler = zeros(size(t_euler));
P_euler(1) = P0; % Initial condition

% Euler Iteration
for i = 2:length(t_euler)
    dPdt = r * P_euler(i-1) * (1 - P_euler(i-1) / K);
    P_euler(i) = P_euler(i-1) + h * dPdt;
end

%% Using ode45
% Define the function for ode45
logistic_growth = @(t, P) r * P * (1 - P / K);
tspan = [t0 tf];

% Solve ODE using ode45
[t_ode, P_ode] = ode45(logistic_growth, tspan, P0);

%% Plot Results
figure;
hold on;
plot(t_euler, P_euler, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(t_ode, P_ode, 'b-', 'LineWidth', 2);
xlabel('Time');
ylabel('Population (P)');
title('Logistic Growth Model: Euler vs. ode45');
legend('Euler Method', 'ode45');
grid on;
hold off;
