clc; clear; close all;

% Parameters
beta = 0.5;  % Infection rate
mu = 0.3;    % Recovery rate
gamma = 0.1; % Immunity loss rate
dt = 0.1;    % Time step for Euler
T = 500;      % Total simulation time
N = 500;     % Number of time steps
t = linspace(0, T, N); % Time vector

% Initial conditions
S0 = 0.9;
I0 = 0.1;
R0 = 0;
y0 = [S0; I0; R0];

%% Euler's Method
S = zeros(1, N);
I = zeros(1, N);
R = zeros(1, N);
S(1) = S0;
I(1) = I0;
R(1) = R0;

for k = 1:N-1
    dS = -beta * S(k) * I(k) + gamma * R(k);
    dI = beta * S(k) * I(k) - mu * I(k);
    dR = mu * I(k) - gamma * R(k);
    
    S(k+1) = S(k) + dt * dS;
    I(k+1) = I(k) + dt * dI;
    R(k+1) = R(k) + dt * dR;
end

%% ODE45 Method
sirs_ode = @(t, y) [-beta * y(1) * y(2) + gamma * y(3); 
                     beta * y(1) * y(2) - mu * y(2);
                     mu * y(2) - gamma * y(3)];
                 
[t_ode, y_ode] = ode45(sirs_ode, [0 T], y0);

%% Plot Results
figure;
hold on;
plot(t, S, 'b--', 'LineWidth', 1.5); % Euler S
plot(t_ode, y_ode(:,1), 'b', 'LineWidth', 2); % ODE45 S

plot(t, I, 'r--', 'LineWidth', 1.5); % Euler I
plot(t_ode, y_ode(:,2), 'r', 'LineWidth', 2); % ODE45 I

plot(t, R, 'g--', 'LineWidth', 1.5); % Euler R
plot(t_ode, y_ode(:,3), 'g', 'LineWidth', 2); % ODE45 R

legend('Euler S', 'ODE45 S', 'Euler I', 'ODE45 I', 'Euler R', 'ODE45 R');
xlabel('Time');
ylabel('Population Fractions');
title('SIRS Model using Euler and ODE45');
grid on;
hold off;
