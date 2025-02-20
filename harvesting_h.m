clc; clear; close all;

%% Parameters
r = 0.2;          % Growth rate
K = 1000;          % Carrying capacity
h = 6;            % Harvesting rate (chosen to cause population drop)
t0 = 0;           % Initial time (not zero)
P0 = 50;          % Initial fish population
tf = 50;          % Final time
h_step = 0.15;     % Step size for Euler's method

%% Euler's Method
t_euler = t0:h_step:tf;  % Time vector
P_euler = zeros(size(t_euler));
P_euler(1) = P0; % Initial condition

% Euler Iteration
for i = 2:length(t_euler)
    dPdt = r * P_euler(i-1) * (1 - P_euler(i-1) / K) - h;
    P_euler(i) = P_euler(i-1) + h_step * dPdt;
    
    % Ensure population does not go negative
    if P_euler(i) < 0
        P_euler(i) = 0;
        break;
    end
end

%% Using ode45
% Define the function for ode45
fish_harvesting = @(t, P) r * P * (1 - P / K) - h;
tspan = [t0 tf];

% Solve ODE using ode45
[t_ode, P_ode] = ode45(fish_harvesting, tspan, P0);

%% Plot Results
figure;
hold on;
plot(t_euler, P_euler, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(t_ode, P_ode, 'b-', 'LineWidth', 2);
xlabel('Time');
ylabel('Fish Population');
title('Fish Harvesting Model: Euler vs. ode45');
legend('Euler Method', 'ode45');
grid on;
hold off;
