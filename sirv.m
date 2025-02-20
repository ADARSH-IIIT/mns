clc; clear; close all;

%% Parameters
beta = 0.3;  % Transmission rate
gamma = 0.1; % Recovery rate
nu = 0.05;   % Vaccination rate
N = 1000;    % Total population
S0 = 990;    % Initial susceptible individuals
I0 = 10;     % Initial infected individuals
R0 = 0;      % Initial recovered individuals
V0 = 0;      % Initial vaccinated individuals

%% Time settings
t0 = 5;     % Initial time (not zero)
t_end = 200; % End time
h = 0.1;    % Step size for Euler

%% ODE45 Solution
sirv_ode = @(t, y) [-beta * y(1) * y(2) / N - nu * y(1); beta * y(1) * y(2) / N - gamma * y(2); gamma * y(2); nu * y(1)];
[t_ode, y_ode] = ode45(sirv_ode, [t0 t_end], [S0; I0; R0; V0]);

%% Euler Method Solution
t_euler = t0:h:t_end;
S_euler = zeros(size(t_euler));
I_euler = zeros(size(t_euler));
R_euler = zeros(size(t_euler));
V_euler = zeros(size(t_euler));
S_euler(1) = S0;
I_euler(1) = I0;
R_euler(1) = R0;
V_euler(1) = V0;

for k = 1:length(t_euler)-1
    dS = -beta * S_euler(k) * I_euler(k) / N - nu * S_euler(k);
    dI = beta * S_euler(k) * I_euler(k) / N - gamma * I_euler(k);
    dR = gamma * I_euler(k);
    dV = nu * S_euler(k);
    
    S_euler(k+1) = S_euler(k) + h * dS;
    I_euler(k+1) = I_euler(k) + h * dI;
    R_euler(k+1) = R_euler(k) + h * dR;
    V_euler(k+1) = V_euler(k) + h * dV;
end

%% Plot Results
figure;
plot(t_ode, y_ode(:,1), 'b', 'LineWidth', 2);
hold on;
plot(t_ode, y_ode(:,2), 'r', 'LineWidth', 2);
plot(t_ode, y_ode(:,3), 'g', 'LineWidth', 2);
plot(t_ode, y_ode(:,4), 'm', 'LineWidth', 2);
plot(t_euler, S_euler, 'b--', 'LineWidth', 2);
plot(t_euler, I_euler, 'r--', 'LineWidth', 2);
plot(t_euler, R_euler, 'g--', 'LineWidth', 2);
plot(t_euler, V_euler, 'm--', 'LineWidth', 2);
hold off;
legend('Susceptible (ODE45)', 'Infected (ODE45)', 'Recovered (ODE45)', 'Vaccinated (ODE45)', 'Susceptible (Euler)', 'Infected (Euler)', 'Recovered (Euler)', 'Vaccinated (Euler)');
xlabel('Time');
ylabel('Population');
title('Susceptible-Infected-Recovered-Vaccinated (SIRV) Model');
grid on;