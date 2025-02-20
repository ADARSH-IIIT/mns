clc; clear; close all;

%% Parameters
r = 0.5;          % Growth rate
K = 100;          % Carrying capacity
h = 0.1;          % Proportional harvesting factor
t0 = 2;           % Initial time (not zero)
P0 = 40;          % Initial fish population
tf = 15 ;          % Final time
h_step = 0.5;     % Step size for Euler's method

%% Euler's Method
t_euler = t0:h_step:tf;  % Time vector
P_euler = zeros(size(t_euler));
P_euler(1) = P0; % Initial condition

% Euler Iteration
for i = 2:length(t_euler)
    dPdt = r * P_euler(i-1) * (1 - P_euler(i-1) / K) - h * P_euler(i-1);
    P_euler(i) = P_euler(i-1) + h_step * dPdt;
    
    % Ensure population does not go negative
    if P_euler(i) < 0
        P_euler(i) = 0;
        t_euler = t_euler(1:i);  % Trim time vector to match
        break;
    end
end

%% Using ode45
% Define the function for ode45
fish_harvesting = @(t, P) r * P * (1 - P / K) - h * P;
tspan = [t0 tf];

% Solve ODE using ode45
[t_ode, P_ode] = ode45(fish_harvesting, tspan, P0);

% Ensure population never goes negative
P_ode(P_ode < 0) = 0;

%% Plot Results
figure;
hold on;
plot(t_euler, P_euler, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(t_ode, P_ode, 'b-', 'LineWidth', 2);
xlabel('Time');
ylabel('Fish Population');
title('Fish Harvesting Model: Proportional Harvesting');
legend('Euler Method', 'ode45');
grid on;
hold off;
