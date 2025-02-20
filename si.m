clc; clear; close all;

%% Parameters
beta = 0.3;  % Transmission rate
N = 1000;    % Total population
S0 = 990;    % Initial susceptible individuals
I0 = 10;     % Initial infected individuals

%% Time settings
t0 = 5;     % Initial time (not zero)
t_end = 500; % End time
h = 0.1;    % Step size for Euler

%% ODE45 Solution
sir_ode = @(t, y) [-beta * y(1) * y(2)/N ; beta * y(1) * y(2)/N ];
[t_ode, y_ode] = ode45(sir_ode, [t0 t_end], [S0; I0]);

%% Euler Method Solution
t_euler = t0:h:t_end;
S_euler = zeros(size(t_euler));
I_euler = zeros(size(t_euler));
S_euler(1) = S0;
I_euler(1) = I0;

for k = 1:length(t_euler)-1
    dS = -beta * S_euler(k) * I_euler(k) / N;
    dI = beta * S_euler(k) * I_euler(k) / N;
    
    S_euler(k+1) = S_euler(k) + h * dS;
    I_euler(k+1) = I_euler(k) + h * dI;
end

%% Plot Results
figure;
plot(t_ode, y_ode(:,1), 'b', 'LineWidth', 2);
hold on;
plot(t_ode, y_ode(:,2), 'r', 'LineWidth', 2);
plot(t_euler, S_euler, 'b--', 'LineWidth', 2);
plot(t_euler, I_euler, 'r--', 'LineWidth', 2);
hold off;
legend('Susceptible (ODE45)', 'Infected (ODE45)', 'Susceptible (Euler)', 'Infected (Euler)');
xlabel('Time');
ylabel('Population');
title('Susceptible-Infected Model');
grid on;