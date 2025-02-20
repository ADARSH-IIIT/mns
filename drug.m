% This MATLAB script simulates a Two-Compartment Pharmacokinetic Model 
% using the ODE45 solver. The model describes how a drug moves between 
% the bloodstream and tissues and how it is eliminated over time.

% --------------------------
% 1. Simulation Overview
% --------------------------
% We are simulating the concentration of a drug in two compartments:
% - Bloodstream (Central Compartment): Where the drug is administered.
% - Tissues (Peripheral Compartment): Where the drug is absorbed and later returned to the blood.
%
% The drug is introduced via IV infusion at a constant rate.
% The drug undergoes:
% - Absorption from the bloodstream to tissues at rate 'ka'.
% - Return from tissues to the bloodstream at rate 'kr'.
% - Elimination from the bloodstream at rate 'ke'.
%
% To analyze drug concentration over time, we solve this system using ODE45, 
% a built-in MATLAB solver for ordinary differential equations.

% --------------------------
% 2. Differential Equations
% --------------------------
% The drug concentrations in the bloodstream (Cb) and tissues (Ct) are modeled 
% using the following first-order differential equations:

% dCb/dt = (R/Vb) - (ka * Cb) + (kr * Ct) - (ke * Cb)
% dCt/dt = (ka * Cb) - (kr * Ct)

% --------------------------
% 3. Parameter Definitions
% --------------------------
% Cb  -> Drug concentration in the bloodstream (mg/L)
% Ct  -> Drug concentration in tissues (mg/L)
% R   -> IV infusion rate (mg/hour)
% Vb  -> Volume of blood (liters)
% Vt  -> Volume of tissues (liters)
% ka  -> Rate of absorption from blood to tissues (/hour)
% kr  -> Rate of return from tissues to blood (/hour)
% ke  -> Elimination rate from blood (/hour)

% --------------------------
% 4. Breakdown of Equations
% --------------------------

% Bloodstream Equation (dCb/dt):
% - (R / Vb)  -> Drug infusion into the bloodstream.
% - (- ka * Cb)  -> Drug moving to tissues.
% - (+ kr * Ct)  -> Drug returning from tissues.
% - (- ke * Cb)  -> Drug being eliminated from the body.

% Tissues Equation (dCt/dt):
% - (+ ka * Cb)  -> Drug absorbed from blood into tissues.
% - (- kr * Ct)  -> Drug returning to blood from tissues.

% --------------------------
% 5. Solving the Equations
% --------------------------
% Using the ODE45 solver, we integrate these equations over a given 
% time span (e.g., 48 hours) to observe how drug concentrations 
% change in both compartments.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Two-Compartment Pharmacokinetic Model
% Author: [Your Name]
% Date:   [Date]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% ----------------------------
% 1) Define model parameters
% ----------------------------
R   = 50;       % IV infusion rate (mg/hour)
Vb  = 5;        % Volume of bloodstream (liters)
Vt  = 10;       % Volume of tissues (liters)
ka  = 0.2;      % Absorption rate from blood to tissues (1/hour)
kr  = 0.1;      % Return rate from tissues to blood (1/hour)
ke  = 0.05;     % Elimination rate from bloodstream (1/hour)

% Initial conditions (assume no drug in body initially)
Cb0 = 0;        % Initial drug concentration in bloodstream (mg/L)
Ct0 = 0;        % Initial drug concentration in tissues (mg/L)

% Time span for simulation
tSpan = [0 48]; % Simulate for 48 hours

% ----------------------------
% 2) Pack initial conditions
% ----------------------------
y0 = [Cb0; Ct0];

% ----------------------------
% 3) Solve the system using ODE45
% ----------------------------
[tSol, ySol] = ode45(@(t,y) drugModel(t, y, R, Vb, Vt, ka, kr, ke), tSpan, y0);

% Extract results
Cb = ySol(:,1);  % Drug concentration in bloodstream
Ct = ySol(:,2);  % Drug concentration in tissues

% ----------------------------
% 4) Plot results
% ----------------------------
figure;
plot(tSol, Cb, 'b', 'LineWidth', 2); hold on;
plot(tSol, Ct, 'r', 'LineWidth', 2);
xlabel('Time (hours)');
ylabel('Drug Concentration (mg/L)');
legend('Bloodstream (C_b)','Tissues (C_t)','Location','best');
title('Pharmacokinetics: Drug Concentration Over Time');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Supporting function: Two-Compartment Pharmacokinetic Model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dYdt = drugModel(t, Y, R, Vb, Vt, ka, kr, ke)
    % Y(1) = Cb (Bloodstream concentration), Y(2) = Ct (Tissue concentration)
    Cb = Y(1);
    Ct = Y(2);

    % Equations:
    dCb_dt = (R / Vb) - ka * Cb + kr * Ct - ke * Cb;
    dCt_dt = ka * Cb - kr * Ct;

    % Return a column vector
    dYdt = [dCb_dt; dCt_dt];
end