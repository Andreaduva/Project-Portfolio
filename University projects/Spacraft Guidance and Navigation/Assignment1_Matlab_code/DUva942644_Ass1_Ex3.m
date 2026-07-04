%% ASSIGNMENT 1 EX 3

% initialize kernels and other settings
close all;
clear; clc;
cspice_kclear();
% restricted 2 body model

% change plots appearance 
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% load spice kernels from meta-kernel - NB to be run while in the Assignment01 folder as current
% working directory
cspice_furnsh('..\assignment01_win.tm');  % Windows version
% cspice_furnsh('./assignment01_mac.tm');  % MacOs version
kernel_number = cspice_ktotal('ALL');
fprintf('\nNumber of kernels loaded is %d\n',kernel_number);

%% EX 3 POINT 2
% set units for length and time (mass still in kg)
L = 149597870.7;    % [km/AU]
T1 = 86400;          % [s/day]

% retrieve useful parameters expressed in the chosen units [AU,days,kg]
mu_Sun = cspice_bodvrd('SUN','GM',1)/(L^3)*(T1^2);        % Sun's gravitational parameter   [AU^3/day^2]
frame = 'ECLIPJ2000';                                    % Using the Ecliptic inertial frame reference frame
observer = 'SUN';                                        % Coordinates taken with respect to the Sun
t_ref = cspice_str2et('2022-08-03-12:45:20 UTC');        % [s]
t_dep1 = t_ref/T1;      % fixed departure date from Earth [days]
m0 = 1500;                                               % S/C initial mass        [kg]
T_eng = 150e-6/L*(T1^2);                                  % max thrust per engine   [kg AU / day^2]
I_sp = 3000/T1;                                           % engine specific impulse [days]
N1 = 4;                                                  % number of thrusters
g0 = 9.81e-3/L*(T1^2);                                    % gravitational fied at sea level [AU/day^2]
parameters1 = [T_eng; N1; I_sp; g0];
units1 = [t_dep1;L;T1];

% make first integration with guessed values and plot it for graphical inspection
tof_guess1 = 6*31;   % time of flight [days]
t_arr_guess1 = t_dep1 + tof_guess1;     % first guess of 5 months of travel [days]
Earth_state = cspice_spkezr('EARTH',t_dep1*T1,frame,'NONE',observer);
E_pos_t0 = Earth_state(1:3)/L;   % [AU]
E_vel_t0 = Earth_state(4:6)/L*T1; % [AU/day]
x0 = [E_pos_t0; E_vel_t0; m0];

% create struct to be fed as input to defined functions
data1 = struct('mu',mu_Sun,'frame',frame,'observer',observer,'parameters',parameters1,'units',units1);
% lambda guesses:
lambda_r_guess1 = [1;1;1]*1e2;
lambda_v_guess1 = -[1;1;1]*1e5;       % TIME IN DAYS
lambda_m_guess1 = 50; 


y0_guess1 = [x0; lambda_r_guess1; lambda_v_guess1; lambda_m_guess1];
% simple propagation using guess values
opt_ode = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
[tt_guess1,yy_guess1] = ode113(@(t,y) TPBVP_ode(t,y,data1), [0 tof_guess1], y0_guess1, opt_ode);
% plot results
%handles_guess1 = plot_routine(tt_guess1,yy_guess1,data1);

% plot guess trajectory
Mars_state = cspice_spkezr('MARS',t_arr_guess1*T1,data1.frame,'NONE',data1.observer);
M_pos_t0 = Mars_state(1:3)/L;   % [AU]
EM_vel_t0 = Mars_state(4:6)/L*T1; % [AU/day]
figure();
grid on;
hold on;
scatter3(0,0,0,85,[0.9294 0.8039 0.3803],'o','filled');
scatter3(E_pos_t0(1),E_pos_t0(2),E_pos_t0(3),85,[0.467 0.674 0.188],'o','filled');
scatter3(M_pos_t0(1),M_pos_t0(2),M_pos_t0(3),85,[1 0.411764705882353  0.16078431372549],'o','filled');
plot3(yy_guess1(:,1), yy_guess1(:,2), yy_guess1(:,3),'Color',"#0072BD",'LineWidth',2);
legend('Sun','Earth','Mars','trajectory','FontSize',15);
xlabel('$x$ [AU]'); ylabel('$y$ [AU]'); zlabel('$z$ [AU]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15; ax.ZAxis.FontSize = 15;


%call fsolve to solve the shooting problem:
opt_fsolve = optimoptions('fsolve','Algorithm','trust-region','Display','iter-detailed','FunctionTolerance',1e-14,...
                           'OptimalityTolerance',1e-14,'MaxFunctionEvaluations',40000,'MaxIterations',4000,'StepTolerance',1e-9);
lambda0_guess1 = [lambda_r_guess1; lambda_v_guess1; lambda_m_guess1; t_arr_guess1];
F0 = shooting_function(lambda0_guess1,x0,data1);
disp(F0);
disp(norm(F0)^2);
tic;
[lambda0_sol1,fval1,flag1,~] = fsolve( @(lambda0) shooting_function(lambda0,x0,data1),lambda0_guess1,opt_fsolve);
toc;

% elaborate on results of fsolve
tof_sol1 = lambda0_sol1(8) - t_dep1;                      % already in days
tf_sol1 = lambda0_sol1(8);
day_sol1 = cspice_et2utc(lambda0_sol1(8)*T1,'ISOC',3);   % calendar arrival day
y0_sol1 = [x0; lambda0_sol1(1:7)];                       % initial state
lambda0_sol1 = lambda0_sol1(1:7);                        % lambda at t_0
[tt_sol1,yy_sol1] = ode113(@(t,y) TPBVP_ode(t,y,data1), [0 tof_sol1], y0_sol1, opt_ode); % integrate to get trajectory
[u_sol1,alpha_sol1,St_sol1] = control_action(yy_sol1,data1);   % thrust factor and direction angle

% Result plots:
E_pos_tf = cspice_spkpos('EARTH',(t_dep1 + tt_sol1(1))*T1,data1.frame,'NONE',data1.observer)/L;  % [AU]
Mars_state = cspice_spkezr('MARS',(t_dep1 + tt_sol1(end))*T1,data1.frame,'NONE',data1.observer); 
M_pos_tf = Mars_state(1:3)/L;                                                 % [AU]
M_vel_tf = Mars_state(4:6)/L*T1;                                              % [AU/day]

% plot solution trajectory
figure();
grid on;
hold on;
scatter3(0,0,0,85,[0.9294 0.8039 0.3803],'o','filled');
scatter3(E_pos_tf(1),E_pos_tf(2),E_pos_tf(3),85,[0.467 0.674 0.188],'o','filled');
scatter3(M_pos_tf(1),M_pos_tf(2),M_pos_tf(3),85,[1 0.411764705882353  0.16078431372549],'o','filled');
plot3(yy_sol1(:,1), yy_sol1(:,2), yy_sol1(:,3),'Color',"#0072BD",'LineWidth',2);
legend('Sun','Earth','Mars','trajectory','FontSize',15);
xlabel('$x$ [AU]'); ylabel('$y$ [AU]'); zlabel('$z$ [AU]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15; ax.ZAxis.FontSize = 15;

% plot behaviour of the states
% POSITION X
figure(); grid on; hold on;
plot(tt_sol1,yy_sol1(:,1)*L,'LineWidth',2);
scatter(tt_sol1(end),M_pos_tf(1)*L,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$x$ $\left[km\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol1(1) tt_sol1(end)]);

% POSITION Y
figure(); grid on; hold on;
plot(tt_sol1,yy_sol1(:,2)*L,'LineWidth',2);
scatter(tt_sol1(end),M_pos_tf(2)*L,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$y$ $\left[km\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol1(1) tt_sol1(end)]);

% POSITION Z
figure(); grid on; hold on;
plot(tt_sol1,yy_sol1(:,3)*L,'LineWidth',2);
scatter(tt_sol1(end),M_pos_tf(3)*L,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$z$ $\left[km\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol1(1) tt_sol1(end)]);

% VELOCITY VX
figure(); grid on; hold on;
plot(tt_sol1,yy_sol1(:,4)*L/T1,'LineWidth',2);
scatter(tt_sol1(end),M_vel_tf(1)*L/T1,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$v_x$ $\left[\frac{km}{s}\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol1(1) tt_sol1(end)]);

% VELOCITY VY
figure(); grid on; hold on;
plot(tt_sol1,yy_sol1(:,5)*L/T1,'LineWidth',2);
scatter(tt_sol1(end),M_vel_tf(2)*L/T1,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$v_y$ $\left[\frac{km}{s}\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol1(1) tt_sol1(end)]);

% VELOCITY VZ
figure(); grid on; hold on;
plot(tt_sol1,yy_sol1(:,6)*L/T1,'LineWidth',2);
scatter(tt_sol1(end),M_vel_tf(3)*L/T1,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$v_x$ $\left[\frac{km}{s}\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol1(1) tt_sol1(end)]);

% LAMBDA M
figure(); grid on; hold on;
plot(tt_sol1,yy_sol1(:,14),'LineWidth',2);
scatter(tt_sol1(end),0,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$\lambda_m$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol1(1) tt_sol1(end)]);

% CONTROL ACTION U
figure(); grid on; hold on;
plot(tt_sol1,u_sol1(:),'LineWidth',2);
xlabel('$t$ [days]'); ylabel('thrust factor [-]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol1(1) tt_sol1(end)]); pbaspect([3 1 1]); 

% SWITCHING FUNCTION U
figure(); grid on; hold on;
plot(tt_sol1,St_sol1(:),'Color',"#A2142F",'LineWidth',2);
xlabel('$t$ [days]'); ylabel('Switching function [-]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol1(1) tt_sol1(end)]);
pbaspect([3 1 1]); 

% MASS
figure(); grid on; hold on;
plot(tt_sol1,yy_sol1(:,7),'Color',"#D95319",'LineWidth',2);
xlabel('$t$ [days]'); ylabel('$m$ $\left[kg\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol1(1) tt_sol1(end)]);
ylim([0 m0]); pbaspect([3 1 1]);


% TRAJECTORY WITH THRUST VECTOR
index_samples = ceil(linspace(1,length(tt_sol1),150));
origin_vec = yy_sol1(index_samples,1:3);
alpha_vec = alpha_sol1(index_samples,:);
figure();
grid on;
hold on;
scatter3(0,0,0,85,[0.9294 0.8039 0.3803],'o','filled');
scatter3(E_pos_tf(1),E_pos_tf(2),E_pos_tf(3),85,[0.467 0.674 0.188],'o','filled');
scatter3(M_pos_tf(1),M_pos_tf(2),M_pos_tf(3),85,[1 0.411764705882353  0.16078431372549],'o','filled');
plot3(yy_sol1(:,1), yy_sol1(:,2), yy_sol1(:,3),'Color',"#0072BD",'LineWidth',2);
quiver3(origin_vec(:,1),origin_vec(:,2),origin_vec(:,3), ...
    0.2*alpha_vec(:,1),0.2*alpha_vec(:,2),0.02*alpha_vec(:,3),'off','Color',"#D95319");
legend('Sun','Earth','Mars','trajectory','thrust direction','FontSize',15,'Location','southwest');
xlabel('$x$ [AU]'); ylabel('$y$ [AU]'); zlabel('$z$ [AU]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15; ax.ZAxis.FontSize = 15;

%% EX 3 POINT 3
% modify number of thrusters and units:
T2 = cspice_spd()*365;  % [s/year]
m2 = 1500;              % mass unit

mu_Sun2 = cspice_bodvrd('SUN','GM',1)/(L^3)*(T2^2);        % Sun's gravitational parameter   [AU^3/day^2]
t_dep2 = t_ref/T2;
m0 = 1500/m2;                                               % S/C initial mass        [-]
T_eng2 = 150e-6/L*(T2^2)/m2;                                  % max thrust per engine   [- AU / year^2]
I_sp2 = 3000/T2;                                           % engine specific impulse [year]
N2 = 3;                                                   % number of thrusters
g02 = 9.81e-3/L*(T2^2);                                    % gravitational fied at sea level [AU/year^2]
parameters2 = [T_eng2; N2; I_sp2; g02];
units2 = [t_dep2; L; T2];
data2 = struct('mu',mu_Sun2,'frame',frame,'observer',observer,'parameters',parameters2,'units',units2);

% make first integration with guessed values and plot it for graphical inspection
tof_guess2 = 2;   % time of flight [days]
t_arr_guess2 = t_dep2 + tof_guess2;     % first guess of 5 months of travel [days]
Earth_state = cspice_spkezr('EARTH',t_dep2*T2,frame,'NONE',observer);
E_pos_t0 = Earth_state(1:3)/L;   % [AU]
E_vel_t0 = Earth_state(4:6)/L*T2; % [AU/day]
x0 = [E_pos_t0; E_vel_t0; m0];

% lambda guesses:
lambda_r_guess2 = [1;1;10];
lambda_v_guess2 = [1;1;10];       % TIME IN DAYS
lambda_m_guess2 = 1; 


y0_guess2 = [x0; lambda_r_guess2; lambda_v_guess2; lambda_m_guess2];
% simple propagation using guess values
opt_ode = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
[tt_guess2,yy_guess2] = ode113(@(t,y) TPBVP_ode(t,y,data2), [0 tof_guess2], y0_guess2, opt_ode);

% call fsolve 
lambda0_guess2 = [lambda_r_guess2; lambda_v_guess2; lambda_m_guess2; t_arr_guess2];
opt_fsolve = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display', ...
         'iter-detailed','FiniteDifferenceType','central','FunctionTolerance',1e-14,...
                          'OptimalityTolerance',1e-10,'MaxFunctionEvaluations',90000, ...
                          'MaxIterations',50000,'StepTolerance',1e-18);
tic;
[lambda0_sol2,fval2,flag2,~] = fsolve( @(lambda0) shooting_function(lambda0,x0,data2),lambda0_guess2,opt_fsolve);
toc;

% elaborate on results of fsolve
tof_sol2 = lambda0_sol2(8) - t_dep2;                      % already in days
tf_sol2 = lambda0_sol2(8);
day_sol2 = cspice_et2utc(lambda0_sol2(8)*T2,'ISOC',3);   % calendar arrival day
y0_sol2 = [x0; lambda0_sol2(1:7)];                       % initial state
lambda0_sol2 = lambda0_sol2(1:7);                        % lambda at t_0
[tt_sol2,yy_sol2] = ode113(@(t,y) TPBVP_ode(t,y,data2), [0 tof_sol2], y0_sol2, opt_ode); % integrate to get trajectory
[u_sol2,alpha_sol2,St_sol2] = control_action(yy_sol2,data2);   % thrust factor and direction angle

% Result plots:
E_pos_tf = cspice_spkpos('EARTH',(t_dep2 + tt_sol2(1))*T2,data2.frame,'NONE',data2.observer)/L;  % [AU]
Mars_state = cspice_spkezr('MARS',(t_dep2 + tt_sol2(end))*T2,data2.frame,'NONE',data2.observer); 
M_pos_tf = Mars_state(1:3)/L;                                                 % [AU]
M_vel_tf = Mars_state(4:6)/L*T2;                                              % [AU/day]

% plot solution trajectory
figure();
grid on;
hold on;
scatter3(0,0,0,85,[0.9294 0.8039 0.3803],'o','filled');
scatter3(E_pos_tf(1),E_pos_tf(2),E_pos_tf(3),85,[0.467 0.674 0.188],'o','filled');
scatter3(M_pos_tf(1),M_pos_tf(2),M_pos_tf(3),85,[1 0.411764705882353  0.16078431372549],'o','filled');
plot3(yy_sol2(:,1), yy_sol2(:,2), yy_sol2(:,3),'Color',"#0072BD",'LineWidth',2);
legend('Sun','Earth','Mars','trajectory','FontSize',15);
xlabel('$x$ [AU]'); ylabel('$y$ [AU]'); zlabel('$z$ [AU]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15; ax.ZAxis.FontSize = 15;

% plot behaviour of the states
% POSITION X
figure(); grid on; hold on;
plot(tt_sol2*T2/T1,yy_sol2(:,1)*L,'LineWidth',2);
scatter(tt_sol2(end)*T2/T1,M_pos_tf(1)*L,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$x$ $\left[km\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol2(1) tt_sol2(end)]);

% POSITION Y
figure(); grid on; hold on;
plot(tt_sol2*T2/T1,yy_sol2(:,2)*L,'LineWidth',2);
scatter(tt_sol2(end)*T2/T1,M_pos_tf(2)*L,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$y$ $\left[km\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol2(1) tt_sol2(end)]);

% POSITION Z
figure(); grid on; hold on;
plot(tt_sol2*T2/T1,yy_sol2(:,3)*L,'LineWidth',2);
scatter(tt_sol2(end)*T2/T1,M_pos_tf(3)*L,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$z$ $\left[km\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol2(1)*T2/T1 tt_sol2(end)*T2/T1]);

% VELOCITY VX
figure(); grid on; hold on;
plot(tt_sol2*T2/T1,yy_sol2(:,4)*L/T2,'LineWidth',2);
scatter(tt_sol2(end)*T2/T1,M_vel_tf(1)*L/T2,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$v_x$ $\left[\frac{km}{s}\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol2(1)*T2/T1 tt_sol2(end)*T2/T1]);

% VELOCITY VY
figure(); grid on; hold on;
plot(tt_sol2*T2/T1,yy_sol2(:,5)*L/T2,'LineWidth',2);
scatter(tt_sol2(end)*T2/T1,M_vel_tf(2)*L/T2,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$v_y$ $\left[\frac{km}{s}\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol2(1)*T2/T1 tt_sol2(end)*T2/T1]);

% VELOCITY VZ
figure(); grid on; hold on;
plot(tt_sol2*T2/T1,yy_sol2(:,6)*L/T2,'LineWidth',2);
scatter(tt_sol2(end)*T2/T1,M_vel_tf(3)*L/T2,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$v_x$ $\left[\frac{km}{s}\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol2(1)*T2/T1 tt_sol2(end)*T2/T1]);

% LAMBDA M
figure(); grid on; hold on;
plot(tt_sol2*T2/T1,yy_sol2(:,14),'LineWidth',2);
scatter(tt_sol2(end)*T2/T1,0,65,'o','filled');
xlabel('$t$ [days]'); ylabel('$\lambda_m$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol2(1)*T2/T1 tt_sol2(end)*T2/T1]);

% CONTROL ACTION U
figure(); grid on; hold on;
plot(tt_sol2*T2/T1,u_sol2(:),'LineWidth',2);
xlabel('$t$ [days]'); ylabel('thrust factor [-]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol2(1)*T2/T1 tt_sol2(end)*T2/T1]); pbaspect([3 1 1]); 

% SWITCHING FUNCTION U
figure(); grid on; hold on;
plot(tt_sol2*T2/T1,St_sol2(:),'Color',"#A2142F",'LineWidth',2);
xlabel('$t$ [days]'); ylabel('Switching function [-]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol2(1)*T2/T1 tt_sol2(end)*T2/T1]);
pbaspect([3 1 1]); 

% MASS
figure(); grid on; hold on;
plot(tt_sol2*T2/T1,yy_sol2(:,7)*m2,'Color',"#D95319",'LineWidth',2);
xlabel('$t$ [days]'); ylabel('$m$ $\left[kg\right]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlim([tt_sol2(1)*T2/T1 tt_sol2(end)*T2/T1]);
ylim([0 m0*m2]); pbaspect([3 1 1]); 


% TRAJECTORY WITH THRUST VECTOR
index_samples = ceil(linspace(1,length(tt_sol2),500));
origin_vec = yy_sol2(index_samples,1:3);
alpha_vec = alpha_sol2(index_samples,:);
figure();
grid on;
hold on;
scatter3(0,0,0,85,[0.9294 0.8039 0.3803],'o','filled');
scatter3(E_pos_tf(1),E_pos_tf(2),E_pos_tf(3),85,[0.467 0.674 0.188],'o','filled');
scatter3(M_pos_tf(1),M_pos_tf(2),M_pos_tf(3),85,[1 0.411764705882353  0.16078431372549],'o','filled');
plot3(yy_sol2(:,1), yy_sol2(:,2), yy_sol2(:,3),'Color',"#0072BD",'LineWidth',2);
quiver3(origin_vec(:,1),origin_vec(:,2),origin_vec(:,3), ...
    0.2*alpha_vec(:,1),0.2*alpha_vec(:,2),0.02*alpha_vec(:,3),'off','Color',"#D95319");
legend('Sun','Earth','Mars','trajectory','thrust direction','FontSize',15,'Location','southwest');
xlabel('$x$ [AU]'); ylabel('$y$ [AU]'); zlabel('$z$ [AU]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15; ax.ZAxis.FontSize = 15;


% clear kernel pool
cspice_kclear();
%% POINT 2 FUNCTIONS


function Dy = TPBVP_ode(~,y,data)
% Function expressing the differential equations ruling the two point boundary value problem we have
% to solve in the optimization strategy. Consists of 14 differential equations related to the spacecraft
% state and equivalent Lagrange multipliers.
% S/C modelled under assumptions of low-thrust and 2 body dynamics.
%
% INPUT: 
% t  = current time, not used since system is autonomous. 
% y  = 14x1 array of the state, defined as (r,v,m,lambda).
% data = struct containing known values. Made by following fields:
%        mu    = Sun's gravitational parameter.                                  [AU^3/day^2]
%        frame = frame in which coordinates are represented.
%        obs   = observerving body.  
%        parameters = vector of fixed quantities. Defined as:
%              T_eng = maximum thrust per engine
%              N     = number of thrusters
%              I_sp  = thruster specific impulse
%              g0    = Earth's gravitational field at sea level
%        units = vector containing relation between default units and chosen ones:
%              t_dep = departure date in chosen units. Needed to call cspice.    [days]
%              L     = length units conversion factor.                           [km/AU]
%              T     = time units conversion factor.                             [sec/day]
% OUTPUT:
% Dy = current derivative of y, still 14x1 vector.
%

% check size of the parameters vector:
if length(data.parameters) ~= 4
    error('\nThe parameters vector passed to the TPBVP ode must have 4 elements\n');
end

% extract various parameters:
r_vect   = y(1:3);             % position vector
v_vect   = y(4:6);             % velocity vector
m        = y(7);               % mass
lambda_r = y(8:10);            % Lag. mult. for r
lambda_v = y(11:13);           % Lag. mult. for v
lambda_m = y(14);              % Lag. mult. for m

mu     = data.mu;               % gravitational parameter                   [AU^3/day^2]
T_eng  = data.parameters(1);    % maximum thrust per engine                 [kg AU/day]
N      = data.parameters(2);    % number of thrusters
I_sp   = data.parameters(3);    % thruster specific impulse                 [day]
g0     = data.parameters(4);    % Earth's gravitational field at sea level  [AU/day^2]
T_max  = N*T_eng;               % total available thrust                    [kg AU/day]
r_norm = norm(r_vect);          % norm of position value
% compute value of switching function:
St = - norm(lambda_v)*I_sp*g0/m - lambda_m;
if St < 0
    u = 1;  % max throttle
elseif St > 0
    u = 0;  % thrust off
elseif St == 0
    u = 0.5; % actually undefined but useful to detect if St = 0 appears
end

% compute current derivatives:
Dy = zeros(14,1);
Dy(1:3)   = v_vect;                                                                % r_dot
Dy(4:6)   = - mu/(r_norm^3)*r_vect - lambda_v/norm(lambda_v)*T_max*u/m;                % v_dot
Dy(7)     = -u*T_max/(I_sp*g0);                                                        % m_dot
Dy(8:10)  = mu/(r_norm^3)*lambda_v - 3*mu/(r_norm^5)*dot(lambda_v,r_vect)*r_vect;  % lambda_r dot
Dy(11:13) = - lambda_r;                                                            % lambda_v dot
Dy(14)    = -norm(lambda_v)*u*T_max/(m^2);                                             % lambda_m dot
end




function F = shooting_function(y,x0,data)
% function that provides the expression that has to be set equal to zero in order to solve our
% TPBVP. It consists of 8 nonlinear equations in 8 unknows to be solver by numerical methods like
% the Newton's one.
%
% INPUT:
% y    = 8x1 vector of unknowns made by (lambda_r(t0),lambda_v(t0),lambda_m(t0),t_final)
% x0   = 7x1 vector of known initial conditions made by (r(t0),v(t0),m(t0))
% data = struct containing known values. Made by following fields:
%        mu         = Sun's gravitational parameter.                             [AU^3/day^2]
%        frame      = frame in which coordinates are represented.
%        observer   = observerving body.  
%        parameters = vector of fixed quantities. Defined as:
%              T_eng = maximum thrust per engine
%              N     = number of thrusters
%              I_sp  = thruster specific impulse
%              g0    = Earth's gravitational field at sea level
%        units = vector containing relation between default units and chosen ones:
%              t_dep = departure date in chosen units. Needed to call cspice.    [days]
%              L     = length units conversion factor.                           [km/AU]
%              T     = time units conversion factor.                             [sec/day]
%
% OUTPUT:
% F = vector of 8 algebraic nonlinear equations. They depends on flow of odes.

% check the data input has all the fields:
if numel(fieldnames(data)) ~= 5
    error('\nThe data struct for the shooting function must have 5 fields.\n')
end
% extract input quantities:
mu    = data.mu;               % gravitational parameter                  [AU^3/day^2]
frame = data.frame;            % reference frame
obs   = data.observer;         % observing body/origin
T_eng = data.parameters(1);    % maximum thrust per engine                [kg AU/day]
N     = data.parameters(2);    % number of thrusters      
I_sp  = data.parameters(3);    % thruster specific impulse                [day]
g0    = data.parameters(4);    % Earth's gravitational field at sea level [AU/day^2]
T_max    = N*T_eng;            % total available thrust                   [kg AU/day]
t_dep = data.units(1);         %  departure ephemeris time                [days]  
L     = data.units(2);         %                                          [km/AU]
T     = data.units(3);         %                                          [sec/day]

t_f_guess      = y(8);
tof            = t_f_guess - t_dep;
% tof = y(8);
% t_f_guess = t_dep + tof;
% propagate the equations of motion using the current values:
opt = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
state_0 = [x0;y(1:7)];
[~,propagated_state] = ode113(@(t,y) TPBVP_ode(t,y,data), [0 tof], state_0, opt);
% retrieve final states:
r_f = propagated_state(end,1:3)';
v_f = propagated_state(end,4:6)';
m_f = propagated_state(end,7);
lambda_r_f = propagated_state(end,8:10)';
lambda_v_f = propagated_state(end,11:13)';
lambda_m_f = propagated_state(end,14)';
% computeswitching function and thrust factor at final time:
St_f = -norm(lambda_v_f)*(I_sp*g0)/m_f - lambda_m_f;       % switching function at final time
if St_f > 0
    u_f = 0;                                              % thrust factor at final time
elseif St_f < 0
    u_f = 1;
elseif St_f == 0
    u_f = 0.5;
end
% Hamiltonian at final time:
H_f = 1 + dot(lambda_r_f,v_f) + dot(lambda_v_f,-mu*r_f/norm(r_f)) + St_f*u_f*T_max/(I_sp*g0);
% Mars position, velocity and acceleration at final time:
Mars_state_f = cspice_spkezr('MARS',t_f_guess*T,frame,'NONE',obs);
Mars_pos_f = Mars_state_f(1:3)/L;                       % Mars position      [AU]
Mars_vel_f = Mars_state_f(4:6)/L*T;                     % Mars velocity      [AU/day]                                     
Mars_acc_f = -mu/(norm(Mars_pos_f)^3)*Mars_pos_f;       % Mars acceleration  [AU/day^2]

F = zeros(8,1);
F(1:3) = r_f - Mars_pos_f;
F(4:6) = v_f - Mars_vel_f;
F(7)   = lambda_m_f;
F(8)   = H_f - dot(lambda_r_f,Mars_vel_f) - dot(lambda_v_f,Mars_acc_f);
end



function [u,alpha,St_vec] = control_action(y,data)
% function that retrieves the value of the control actions u and alpha at each solution time.
%
% INPUT:
% y  = Nx14 matrix of the state solution in time, defined as (r,v,m,lambda).
%           where N is the number of points in time at which the solution has been computed
% data = struct containing known values. Made by following fields:
%        mu    = Sun's gravitational parameter.                                  [AU^3/day^2]
%        frame = frame in which coordinates are represented.
%        obs   = observerving body.  
%        parameters = vector of fixed quantities. Defined as:
%              T_eng = maximum thrust per engine
%              N     = number of thrusters
%              I_sp  = thruster specific impulse
%              g0    = Earth's gravitational field at sea level
%        units = vector containing relation between default units and chosen ones:
%              t_dep = departure date in chosen units. Needed to call cspice.    [days]
%              L     = length units conversion factor.                           [km/AU]
%              T     = time units conversion factor.                             [sec/day]
% OUTPUT:
% u     = Nx1 array for the throttle factor at each point in time.        [-]
% alpha = Nx3 matrix for the thrust pointing angle at each point in time. [rad]

% check size of the parameters vector:
if length(data.parameters) ~= 4
    error('\nThe parameters vector passed to the TPBVP ode must have 4 elements\n');
end

% extract useful parameters
I_sp  = data.parameters(3);         % thruster specific impulse
g0    = data.parameters(4);         % Earth's gravitational field at sea level
% initialize output variables
[rows,~] = size(y);
St_vec = zeros(rows,1);
u = zeros(rows,1);
alpha = zeros(rows,3);


% start a for loop
for i = 1:rows
    m      = y(i,7);                 % mass
    lambda_v = y(i,11:13);           % Lag. mult. for v
    lambda_m = y(i,14);              % Lag. mult. for m
    % compute value of switching function:
    St = -norm(lambda_v)*I_sp*g0/m - lambda_m;
    if St > 0
        u(i) = 0;                             % thrust throttle
    elseif St < 0
        u(i) = 1;
    else
        u(i) = 5;
    end
    St_vec(i) = St;
    alpha(i,:) = - lambda_v/(norm(lambda_v)); % thrust angle
end
end


function handles = plot_routine(tt,yy,data)
% function that given the integration results plots relevant graphs and
% gives back their handles.
%
% INPUT:
% tt    = Nx1 vector containing times at which solution is computed.          [days]
% yy    = Nx14 matrix containing the augmented state values at each time.
% data = struct containing known values. Made by following fields:
%        mu    = Sun's gravitational parameter.                                  [AU^3/day^2]
%        frame = frame in which coordinates are represented.
%        obs   = observerving body.  
%        parameters = vector of fixed quantities. Defined as:
%              T_eng = maximum thrust per engine
%              N     = number of thrusters
%              I_sp  = thruster specific impulse
%              g0    = Earth's gravitational field at sea level
%              units = vector containing relation between default units and chosen ones:
%                      t_dep = departure date in chosen units. Needed to call cspice.    [days]
%                      L     = length units conversion factor.                           [km/AU]
%                      T     = time units conversion factor.                             [sec/day]

% OUTPUT:
% plots   = 3 figures: one for the spacecraft trajectory, another one for 
%                    time evolution of s/c states and the last for the evolution
%                    of Lagrange Multipliers.
% handles = 15x1 vector containing handles to each resulting plot.

handles = zeros(15,1);
% check size of the units vector:
if length(data.units) ~= 3
    error('\nThe parameters vector passed to the plotting routine must have 3 elements\n');
end
frame = data.frame;           % reference frame
obs   = data.observer;        % observing body/origin
t_dep = data.units(1);        %  departure ephemeris time [days]  
L     = data.units(2);        %  [km/AU]
T     = data.units(3);        %  [sec/day]

Earth_pos = cspice_spkpos('EARTH',(t_dep + tt(1))*T,frame,'NONE',obs)/L;      % [AU]
Mars_state = cspice_spkezr('MARS',(t_dep + tt(end))*T,frame,'NONE',obs); 
Mars_pos = Mars_state(1:3)/L;                                                 % [AU]
Mars_vel = Mars_state(4:6)/L*T;                                               % [AU/day]

% plot trajectory
figure();
grid on;
hold on;
plot3(0,0,0,'o','MarkerSize',10,'MarkerFaceColor','#FF8800','MarkerEdgeColor','k');
plot3(Earth_pos(1),Earth_pos(2),Earth_pos(3),'og','MarkerSize',10,'MarkerFaceColor','g','MarkerEdgeColor','k');
plot3(Mars_pos(1),Mars_pos(2),Mars_pos(3),'or','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','k');
plot3(yy(:,1), yy(:,2), yy(:,3),'b','LineWidth',2);
title('S/C trajectory','FontSize',18);
legend('Sun','Earth','Mars','S/C','FontSize',15);
handles(1) = gca;          % TO CHANGE OR INSPECT LINES PROPERTIES USE GET(HANDLES(i).CHILDREN(j))


% plot behaviour of the states
figure();
subplot(3,3,1); grid on; hold on;
plot(tt,yy(:,1)*L,'LineWidth',2);
plot(tt(end),Mars_pos(1)*L,'or','MarkerSize',5,'MarkerFaceColor','r');
title('x(t)','FontSize',15); handles(2) = gca;

subplot(3,3,2); grid on; hold on;
plot(tt,yy(:,2)*L,'LineWidth',2);
plot(tt(end),Mars_pos(2)*L,'or','MarkerSize',5,'MarkerFaceColor','r');
title('y(t)','FontSize',15); handles(3) = gca;

subplot(3,3,3); grid on; hold on;
plot(tt,yy(:,3)*L,'LineWidth',2);
plot(tt(end),Mars_pos(3)*L,'or','MarkerSize',5,'MarkerFaceColor','r');
title('z(t)','FontSize',15); handles(4) = gca;

subplot(3,3,4); grid on; hold on;
plot(tt,yy(:,4)*L/T,'LineWidth',2);
plot(tt(end),Mars_vel(1)*L/T,'or','MarkerSize',5,'MarkerFaceColor','r');
title('V_x(t)','FontSize',15); handles(5) = gca;

subplot(3,3,5); grid on; hold on;
plot(tt,yy(:,5)*L/T,'LineWidth',2);
plot(tt(end),Mars_vel(2)*L/T,'or','MarkerSize',5,'MarkerFaceColor','r');
title('V_y(t)','FontSize',15); handles(6) = gca;

subplot(3,3,6); grid on; hold on;
plot(tt,yy(:,6)*L/T,'LineWidth',2);
plot(tt(end),Mars_vel(3)*L/T,'or','MarkerSize',5,'MarkerFaceColor','r');
title('V_z(t)','FontSize',15); handles(7) = gca;

subplot(3,3,7); grid on; hold on;
plot(tt,yy(:,7),'LineWidth',2);
%plot(tt(end),0,'or','MarkerSize',5,'MarkerFaceColor','r');
yline(0,'k','LineWidth',2);
title('m(t)','FontSize',15); handles(8) = gca;
sgtitle('S/C state behaviour','FontSize',15);

% plot behaviour of the Lagrange multipliers
figure();
tiledlayout(3,3); nexttile;
grid on; hold on;
plot(tt,yy(:,8),'LineWidth',2);
title('LambdaR_1(t)','FontSize',15); handles(9) = gca;

nexttile;
grid on; hold on;
plot(tt,yy(:,9),'LineWidth',2);
title('LambdaR_2(t)','FontSize',15); handles(10) = gca;

nexttile;
grid on; hold on;
plot(tt,yy(:,10),'LineWidth',2);
title('LambdaR_3(t)','FontSize',15); handles(11) = gca;

nexttile;
grid on; hold on;
plot(tt,yy(:,11),'LineWidth',2);
title('LambdaV_1(t)','FontSize',15); handles(12) = gca;

nexttile;
grid on; hold on;
plot(tt,yy(:,12),'LineWidth',2);
title('LambdaV_2(t)','FontSize',15); handles(13) = gca;

nexttile;
grid on; hold on;
plot(tt,yy(:,13),'LineWidth',2);
title('LambdaV_3(t)','FontSize',15); handles(14) = gca;

nexttile;
grid on; hold on;
plot(tt,yy(:,14),'LineWidth',2);
plot(tt(end),0,'or','MarkerSize',5,'MarkerFaceColor','r');
title('LambdaM(t)','FontSize',15); handles(15) = gca;
sgtitle('Lagrange Multipliers behaviour','FontSize',15);
end
