
clc; clear; close all;

%% Parameters
r = 0.5;         % Growth rate
t0 = 2;         % Initial time (not 0)
y0 = 5;         % Initial value
tf = 10;        % Final time
h = 0.5;        % Step size for Euler's method

%% Euler's Method
t_euler = t0:h:tf;  % Time vector
y_euler = zeros(size(t_euler));
y_euler(1) = y0; % Initial value

% Euler Iteration
for i = 2:length(t_euler)
    y_euler(i) = y_euler(i-1) + h * r * y_euler(i-1);
end

%% Using ode45
% Define the function for ode45
exp_growth = @(t, y) r * y;
tspan = [t0 tf];

% Solve ODE using ode45
[t_ode, y_ode] = ode45(exp_growth, tspan, y0);

%% Plot Results
figure;
hold on;
plot(t_euler, y_euler, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(t_ode, y_ode, 'b-', 'LineWidth', 2);
xlabel('Time');
ylabel('y');
title('Exponential Growth Model: Euler vs. ode45');
legend('Euler Method', 'ode45');
grid on;
hold off;
