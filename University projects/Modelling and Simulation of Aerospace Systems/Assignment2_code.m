%% Andrea D'Uva 942644
close all;
clear;
clc;
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% EX 1

% data
J1 = 0.2; % kg m
J2 = 0.1; % kg m
T0 = 0.1; % N m
t0 = 0; % s
tf = 10; % s

% equations of motions assuming linear damping, expressed as x_dot = f(t,x);
% x = [ theta_1; theta_1_dot; theta_2; thetha_2_dot]
eq_motion_linear = @(k,b,x) [x(2);
                            k/J1*( x(3) - x(1) );
                            x(4);
                            ( k*( x(1) - x(3) ) - b*x(4) + T0 )/J2];

% assume mild steel as material for the shaft and circular section
d_shaft = 0.025; % m
L_shaft = 4; % m
G_shaft = 75e9; % Pa
J_shaft = pi*(d_shaft/2)^4/32; % m^4
k_guess = G_shaft*J_shaft/L_shaft;    % N rad
% assume an arbitrary value for the damping
b_guess = 4;  

% numerical integration through ode45
x0 = [0;0;0;0]; % initial condition
[t_span, x_solution] = ode45( @(t,x) eq_motion_linear(k_guess,b_guess,x), [t0 tf], x0);
% plot the results in the transient region
figure();
plot(t_span,x_solution(:,1), t_span,x_solution(:,2), t_span,x_solution(:,3), t_span,x_solution(:,4), 'LineWidth',2);
legend('$\theta_1$', '$\dot{\theta_1}$', '$\theta_2$', '$\dot{\theta_2}$', 'FontSize', 15);
xlabel('$t [s]$', 'FontSize', 18); ylabel('$[rad], [\frac{rad}{s}]$', 'FontSize', 18);
xlim([0 1.5]); ylim([0 0.05]); grid on; ax = gca; ax.FontSize = 15;

% import the samples
samples = importdata('samples.txt');
time_samples  = samples.data(:,1);
acc_samples = samples.data(:,2:3);

% define a search-grid to minimize the error between samples and numerical
% results
k_span = linspace(0,10,70);
b_span = linspace(0,10,70);
error_grid = zeros(length(k_span), length(b_span));
tic;
for i = 1:length(k_span)
    for j = 1:length(b_span)
        k = k_span(i);
        b = k_span(j);
        [~,solution] = ode45( @(t,x) eq_motion_linear(k,b,x), time_samples, x0);
        acc_num = state_acceleration(@(x) eq_motion_linear(k,b,x), solution);
        error_grid(i,j) = norm( acc_num - acc_samples,2);
        
    end
end
toc;

% find the minimum value
[min_error, row_index, column_index] = matrix_minimum(error_grid);
k_err_min = k_span(row_index);
b_err_min = b_span(column_index);

% plot the resulting grid
figure();
[K,B] = meshgrid(k_span, b_span);
surf(K,B,error_grid', 'EdgeColor','none');
caxis([0 3]);
c = colorbar;
c.Label.String = 'relative error'; c.Label.FontSize = 15;
hold on; 
plot3(k_err_min,b_err_min,min_error+5, '.k','MarkerSize',20,'LineWidth',15);
text_name = "Minimum:"; text_k = sprintf('$k_{min}$ = %.2f', k_err_min);
text_b = sprintf('$b_{min}$ = %.2f', b_err_min); text_err = sprintf('$err_{min}$ = %.2f', min_error);
str = {text_name, text_k, text_b, text_err};
annotation('textbox', [0.38, 0.15, 0.1, 0.1],'String',str,'interpreter','latex', 'BackgroundColor','w','FontSize',15);
view(2);
xlabel('$k$', 'FontSize' , 15); ylabel('$b$', 'FontSize', 15); ax = gca; ax.FontSize = 15;

% plot acceleration of linear model vs samples
[~, solution_best_linear] = ode45( @(t,x) eq_motion_linear(k_err_min,b_err_min,x), time_samples, x0);
acc_num_linear = state_acceleration(@(x) eq_motion_linear(k_err_min,b_err_min,x), solution_best_linear);

% theta 1
figure();
plot(time_samples, acc_samples(:,1),'o'); hold on; grid on;
plot(time_samples, acc_num_linear(:,1), 'LineWidth', 2);
xlabel('$t \left[ s \right]$', 'FontSize' , 18); ylabel('$\left[ \frac{rad}{s^2} \right]$', 'FontSize', 20);
legend('$ \ddot{\theta}_1^{experimental}$','$ \ddot{\theta}_1^{simulated}$', 'FontSize', 18);
ax = gca; ax.FontSize = 15;

% theta 2
figure();
plot(time_samples, acc_samples(:,2),'o'); hold on; grid on;
plot(time_samples, acc_num_linear(:,2), 'LineWidth', 2);
xlabel('$t \left[ s \right]$', 'FontSize' , 18); ylabel('$\left[ \frac{rad}{s^2} \right]$', 'FontSize', 20);
legend('$ \ddot{\theta}_2^{experimental}$','$ \ddot{\theta}_2^{simulated}$', 'FontSize', 18);
ax = gca; ax.FontSize = 15;



% change to non linear damping

eq_motion_nonlinear = @(k,b,x) [ x(2);
                                k/J1*( x(3) - x(1) ) ;
                                 x(4);
                                (k*( x(1) - x(3) )-b*x(4)*abs(x(4)) + T0 )/J2];

% new search-grid with the non-linear damping
k_span = linspace(0,10,70);
b_span = linspace(0,10,70);
error_grid_nl = zeros(length(k_span), length(b_span));
tic;
for i = 1:length(k_span)
    for j = 1:length(b_span)
        k = k_span(i);
        b = k_span(j);
        [~,solution_nl] = ode45( @(t,x) eq_motion_nonlinear(k,b,x), time_samples, x0);
        acc_num_nl = state_acceleration(@(x) eq_motion_nonlinear(k,b,x), solution_nl);
        error_grid_nl(i,j) = norm( acc_num_nl - acc_samples,2);
        
    end
end
toc;

% find the minimum value
[err_min_nl, row_nl, column_nl] = matrix_minimum(error_grid_nl);
k_err_min_nl = k_span(row_nl);
b_err_min_nl = b_span(column_nl);
% plot the resulting grid
figure();
[K,B] = meshgrid(k_span, b_span);
surf(K,B,error_grid_nl', 'EdgeColor','none');
caxis([0 3]);
c = colorbar;
c.Label.String = 'relative error';
c.Label.FontSize = 15;
hold on;
plot3(k_err_min_nl,b_err_min_nl,err_min_nl+5, '.k','MarkerSize',20,'LineWidth',15);
text_name = "Minimum:"; text_k_nl = sprintf('$k_{min}$ = %.2f', k_err_min_nl);
text_b_nl = sprintf('$b_{min}$ = %.2f', b_err_min_nl); text_err_nl = sprintf('$err_{min}$ = %.2f', err_min_nl);
str_nl = {text_name, text_k_nl, text_b_nl, text_err_nl};
annotation('textbox', [0.38, 0.3, 0.1, 0.1],'String',str_nl,'interpreter','latex', 'BackgroundColor','w','FontSize',15);
view(2);
xlabel('$k$', 'FontSize' , 15); ylabel('$b$', 'FontSize', 15);
ax = gca; ax.FontSize = 15;

% plot acceleration of non-linear model vs samples
[~, solution_best] = ode45( @(t,x) eq_motion_nonlinear(k_err_min_nl,b_err_min_nl,x), time_samples, x0);
acc_num_best = state_acceleration(@(x) eq_motion_nonlinear(k_err_min_nl,b_err_min_nl,x), solution_best);

%theta 1
figure();
plot(time_samples, acc_samples(:,1),'o'); hold on; grid on;
plot(time_samples, acc_num_best(:,1), 'LineWidth',2);
xlabel('$t \left[ s \right]$', 'FontSize' , 18); ylabel('$\left[ \frac{rad}{s^2} \right]$', 'FontSize', 20);
legend('$ \ddot{\theta}_1^{experimental}$','$ \ddot{\theta}_1^{simulated}$', 'FontSize', 18);
ax = gca; ax.FontSize = 15;

% theta 2
figure();
plot(time_samples, acc_samples(:,2),'o'); hold on; grid on;
plot(time_samples, acc_num_best(:,2), 'LineWidth', 2);
xlabel('$t \left[ s \right]$', 'FontSize' , 18); ylabel('$\left[ \frac{rad}{s^2} \right]$', 'FontSize', 20);
legend('$ \ddot{\theta}_2^{experimental}$','$ \ddot{\theta}_2^{simulated}$', 'FontSize', 18);
ax = gca; ax.FontSize = 15;

%% EX 2
% set the data
t0 = 0;   % s
%t1 = 1;   % s
t2 = 1.5; % s
tf = 3;   % s
x0 = 0;   % m
v0 = 0;   % m/s

rho = 890; % kg/m^3 assumed constant
% accumulator
accumulator.p_inf = 2.5e6;   % Pa
accumulator.V_inf = 10e-3;    % m^3
accumulator.p0 = 21e6;       % Pa
accumulator.V0 = accumulator.V_inf*accumulator.p_inf/accumulator.p0; % isothermal expansion p*V = const
accumulator.gamma = 1.2; % [-]

% delivery line
delivery.k_acc = 1.12;                  % [-]
delivery.k_cv  = 2;                     % [-]
delivery.D     = 18e-3;                 %  m
delivery.A     = pi/4*(delivery.D)^2;   % m^2
delivery.L     = 2;                     % m
delivery.f     = 0.032;                 % [-]

% distributor
distributor.k  = 12;     % [-] assumed constant 
distributor.r  = 2.5e-3; %  m 

% actuator
actuator.R_cyl =  25e-3;              % m
actuator.R_rod =  11e-3;              % m
actuator.m     =   2;                 % kg
actuator.x_max =  200e-3;             % m
actuator.load  =  @(x) 1e3 + 120e3*x; % N

% return line
rline.D     = 18e-3;                  %  m
rline.A     = pi/4*(rline.D)^2;       %  m^2
rline.L     = 15;                     %  m
rline.f     = 0.035;                  % [-]
rline.kt    = 1.12;                    % [-]

% tank
tank.p = 0.1e6; % Pa
tank.V0 = 1e-3; % m^3

% initial conditions:
x0 = [0; x0; v0; tank.V0];
parameters = {rho, accumulator, delivery, distributor, actuator, rline, tank};

% solve the problem:

% first we find the time at which the actuator stroke occurs, if it's
% before tf we will have trobules in directly integrating in [t0 tf]
opt = odeset('event', @(t,x) actuator_stroke(t,x,actuator.x_max) );
[t_span_stroke, state_to_stroke, t_stroke, x_stroke, ie] = ode23s(@(t,x) hydraulic_network(t,x,parameters), [t0 tf], x0, opt);

% stroke occurs before tf, so we need to split the time domain and setting
% new appropriate initial conditions to have v_actuator = 0; Tried several
% different checks inside the state function but none delivered this
% result.

x0_stroke = [state_to_stroke(end,1), actuator.x_max, 0, state_to_stroke(end,4)];  % initial conditions are variables at stroke except for v_actuator.
[t_span_post_stroke,state_post_stroke] = ode23s(@(t,x) hydraulic_network(t,x,parameters), [t_stroke tf], x0_stroke);

% now we can recover all parameters before and after the stroke by just
% merging the two results of the numerical integration. 
t_span = [t_span_stroke(1:end-1); t_span_post_stroke(:)];              % full integration time
system_state = [state_to_stroke(1:end-1,:); state_post_stroke(:,:)];   % full state matrix
useful_parameters = zeros(length(t_span),16);
for i = 1:length(t_span)
[~,useful_parameters(i,:)] = hydraulic_network(t_span(i),system_state(i,:),parameters);    % for cycle to recover flows and pressures
end


% plots 
figure(); % pressures on delivery line
plot(t_span, useful_parameters(:,9)/1e6, t_span, useful_parameters(:,13)/1e6, 'LineWidth',2);
legend('$P_{accumulator}$','$P_4$', 'FontSize',15);
xlabel('$t \left[ s \right]$', 'FontSize',15); ylabel('Pressure $\left[ MPa \right]$', 'FontSize',15); grid on;
set(gca, 'XTick', sort([t_stroke, get(gca, 'XTick')]));
ax = gca; xlabels = get(ax, 'XTickLabel');
index = find(strcmp(xlabels,'1.8671'));
xlabels{index} = '$t_{stroke}$';   %needs to exist but make it empty
set(ax, 'XTickLabel', xlabels); ax.FontSize = 15;

figure(); % pressures on return line
plot(t_span, useful_parameters(:,14)/1e6, t_span, useful_parameters(:,15)/1e6, 'LineWidth',2);
legend('$P_5$','$P_6$', 'FontSize',15);
xlabel('$t \left[ s \right]$', 'FontSize',15); ylabel('Pressure $\left[ MPa \right]$', 'FontSize',15); grid on;
set(gca, 'XTick', sort([t_stroke, get(gca, 'XTick')]));
set(gca, 'XTickLabel', xlabels); ax = gca; ax.FontSize = 15;

figure(); % flow rates
plot(t_span, useful_parameters(:,1)*1e3, t_span, useful_parameters(:,6)*1e3, 'LineWidth',2);
legend('$Q_{delivery}$','$Q_{return}$', 'FontSize',15);
xlabel('$t \left[ s \right]$', 'FontSize',15); ylabel('Flow rate $\left[ \frac{l}{s} \right]$', 'FontSize',15); grid on;
set(gca, 'XTick', sort([t_stroke, get(gca, 'XTick')]));
set(gca, 'XTickLabel', xlabels); ax = gca; ax.FontSize = 15;

figure(); % non-dimensional comparison
plot(t_span, useful_parameters(:,6)/max(useful_parameters(:,6)), t_span, useful_parameters(:,14)/max(useful_parameters(:,14)), 'LineWidth',2);
legend('$Q_{return}$','$P_5$', 'FontSize',15);
xlabel('$t \left[ s \right]$', 'FontSize',15); ylabel('Normalized amplitude $\left[-\right]$', 'FontSize',15); grid on;
set(gca, 'XTick', sort([t_stroke, get(gca, 'XTick')]));
set(gca, 'XTickLabel', xlabels); ax = gca; ax.FontSize = 15;

figure(); % volumes of N2 and tank
plot(t_span, (accumulator.V0 + system_state(:,1))*1e3, t_span, system_state(:,4)*1e3, 'LineWidth',2);
legend('$V_{N_2}$','$V_{tank}$', 'FontSize',15);
xlabel('$t \left[ s \right]$', 'FontSize',15); ylabel('Volume $\left[ l \right]$', 'FontSize',15); grid on;
set(gca, 'XTick', sort([t_stroke, get(gca, 'XTick')]));
set(gca, 'XTickLabel', xlabels); ax = gca; ax.FontSize = 15;

figure(); % piston movement
plot(t_span, system_state(:,2)*1e3, 'LineWidth', 2);
xlabel('$t \left[ s \right]$', 'FontSize',15); ylabel('$x_{piston} \left[ mm \right]$', 'FontSize',15); grid on;
set(gca, 'XTick', sort([t_stroke, get(gca, 'XTick')]));
set(gca, 'XTickLabel', xlabels); ax = gca; ax.FontSize = 15;

figure(); % piston velocity
plot(t_span, system_state(:,3)*1e3, 'LineWidth', 2);
xlabel('$t \left[ s \right]$', 'FontSize',15); ylabel('$v_{piston} \left[ \frac{mm}{s} \right]$', 'FontSize',15); grid on;
set(gca, 'XTick', sort([t_stroke, get(gca, 'XTick')]));
set(gca, 'XTickLabel', xlabels); ax = gca; ax.FontSize = 15;

%% EX 3

% define data
R1 = 1000; %  ohm
R2 = 100;  %  ohm
L = 1e-3;  %  H
C = 1e-3;  %  F
f = 5;     %  Hz
t0 = 0;    %  s
tf = 3;    %  s

% equation of motion for the capacitor discharge in form x_dot = A*x
% being x = [ V_c ; V_c_dot]
capacitor_discharge = @(x) [    0                            1;
                          -1/(L*C*( 1 + R2/R1))   -(R2/L + 1/(C*R1))/(1 + R2/R1)]*[ x(1); x(2)];
% initial conditions
V0 = [1 ; 0];
% compute eigenvalues of the state matrix to check their position and plot them                     
lambda = eig( [    0                            1;
           -1/(L*C*( 1 + R2/R1))   -(R2/L + 1/(C*R1))/(1 + R2/R1)] );
figure();
plot(lambda(1),0,'xr',lambda(2),0,'xb','MarkerSize',15,'LineWidth',10); hold on; grid on;
xline(0,'LineWidth',2); yline(0,'LineWidth',2);
xlabel('Re$(\lambda)$','FontSize',16); ylabel('Im$(\lambda)$','FontSize',16);
legend('$\lambda_1$','$\lambda_2$','FontSize',15); ax = gca; ax.FontSize = 15;
                      
% Solve numerically the capacitor discharge using an implicit method
% because the eigenvalues hare quite far from the origin.
% An explicit method would require many more steps to solve this problem
[t_discharge, voltage_discharge] = ode23s( @(t,x) capacitor_discharge(x), [t0 tf], V0 );

% plot the solution
figure();
plot(t_discharge, voltage_discharge(:,1), 'r', 'LineWidth', 2); grid on;
xlabel('$t \left[ s \right]$','FontSize',16); ylabel('$V_c \left[ V \right]$','FontSize',16); ylim([0 1]);
legend('$V_c \left( t \right)$','FontSize',16); ax = gca; ax.FontSize = 15;

% case with a oscillating voltage source

% source function and its derivative
V_t = @(t) sin(2*pi*f*t).*atan(t);
V_t_dot = @(t) 2*pi*f*atan(t).*cos(2*pi*f*t) + 1/(1+t^2).*sin(2*pi*f*t);

% equation of motion for the capacitor voltage in form x_dot = A*x + u
% being x = [ V_c ; V_c_dot]
capacitor_forced = @(t,x) [    0                            1;
                          -1/(L*C*( 1 + R2/R1))   -(R2/L + 1/(C*R1))/(1 + R2/R1)]*[ x(1); x(2)] + [0; ( V_t(t)/(L*C) - V_t_dot(t)/(R1*C)*(1 + R2/R1) )];
% assume homogenous initial conditions on the capacitor
V0_forced = [0, 0];
% numerical integration
[t_forced, voltage_forced] = ode23s( @(t,x) capacitor_forced(t,x), [0 tf], V0_forced);

% plot the forced solution with homogeneous initial conditions
figure();
plot(t_forced, voltage_forced(:,1), 'r', 'LineWidth', 2);
hold on; grid on;
plot( t_forced, V_t(t_forced),'b');
xlabel('$t \left[ s \right]$','FontSize',16); ylabel('$V_c \left[ V \right]$','FontSize',16);
legend('$V_c \left( t \right)$', '$V_0 \left( t \right)$','FontSize',16);
ax = gca; ax.FontSize = 15;

% forced solution with non-homogeneous initial conditions
[t_forced_nh, voltage_forced_nh] = ode23s( @(t,x) capacitor_forced(t,x), [0 tf], [1, 0]);
figure();
plot(t_forced_nh, voltage_forced_nh(:,1), 'r', 'LineWidth', 2); grid on;
xlabel('$t \left[ s \right]$','FontSize',16); ylabel('$V_c \left[ V \right]$','FontSize',16);
legend('$V_c \left( t \right)$','FontSize',16); ax = gca; ax.FontSize = 15;

%% EX4 

% set data:
density = [ 0; 1500; 0; 150; 0];           % kg/m^3
conductivity = [ 30; 100; 30; 0.015; 30];  % W/(m K)
specific_h = [ 0 ; 750; 0; 1000; 0];       % J/(kg K)
layer_l = [7; 50; 7; 50; 7]*1e-3;          % m

t0 = 0;  % s
t1 = 1;  % s
tf = 60; % s
T0 = 293.15;   % K 
Tf = 1273.15;  % K

parameters = {density,conductivity, specific_h,layer_l, T0, Tf};

% simulation with 2 nodes
x0_2nodes = [T0 ; T0];
[t_span_2nodes, T_12] = ode45(@(t,x) simulation_2nodes(t,x,parameters), t0:1:tf, x0_2nodes);
% retrieve T_i
T_i = imposed_temperature(t_span_2nodes, T0, Tf);
% plot temporal evolution
figure();
plot(t_span_2nodes, T_i, '-r', t_span_2nodes, T_12, 'LineWidth',2); grid on;
legend('$T_i(t)$', '$T_{node1}$', '$T_{node2}$', 'FontSize',18);
xlabel('$t \left[ s \right]$', 'FontSize',18); ylabel('Temperature $ \left[ K \right]$', 'FontSize',18);
ax = gca; ax.FontSize = 15;

% simulation with 4 nodes
x0_4nodes = [T0; T0; T0; T0];
[t_span_4nodes, T_1234] = ode45(@(t,x) simulation_4nodes(t,x,parameters), t0:1:tf, x0_4nodes);
% retrieve T_i
T_i4nodes = imposed_temperature(t_span_4nodes, T0, Tf);
% plot temporal evolution
figure();
plot(t_span_4nodes, T_i4nodes, 'r', t_span_4nodes, T_1234, 'LineWidth',2); grid on;
legend('$T_i(t)$', '$T_{node1}$', '$T_{node2}$','$T_{node3}$','$T_{node4}$', 'FontSize',18);
xlabel('$t \left[ s \right]$', 'FontSize',18); ylabel('Temperature $ \left[ K \right]$', 'FontSize',18);
ax = gca; ax.FontSize = 15;

% plot temperature profiles at t = tf for 2 nodes
figure();
distance_1 = layer_l(1);
distance_2 = distance_1 + layer_l(2);
distance_3 = distance_2 + layer_l(3);
distance_4 = distance_3 + layer_l(4);
distance_5 = sum(layer_l);
x_points = [0, distance_1, distance_2, distance_3, distance_4, distance_5];
T_points = [Tf, T_12(end,1), T_12(end,1), T_12(end,2), T_12(end,2), T0];
plot(x_points, T_points, 'LineWidth',2); grid on; xlim([x_points(1) x_points(end)]);
xlabel('$x \left[ mm \right]$', 'FontSize',18); ylabel('Temperature $ \left[ K \right]$', 'FontSize',18);
xticks(x_points); 
xticklabels({'0','$l_1$', '$l_2$', '$l_3$', '$l_4$', '$l_5$'});
ax = gca; ax.FontSize = 15; 


% plot temperature profile at t = tf for 4 nodes
figure();
distance_1 = layer_l(1);
distance_2 = distance_1 + layer_l(2)/2;
distance_3 = distance_2 + layer_l(2)/2;
distance_4 = distance_3 + layer_l(3);
distance_5 = distance_4 + layer_l(4)/2;
distance_6 = distance_5 + layer_l(4)/2;
distance_7 = sum(layer_l);
x_points_4nodes = [0, distance_1, distance_2, distance_2, distance_3, distance_4, distance_5, distance_5, distance_6, distance_7];
T_points_4nodes = [Tf, T_1234(end,1), T_1234(end,1), T_1234(end,2), T_1234(end,2), T_1234(end,3),T_1234(end,3),T_1234(end,4),T_1234(end,4), T0];
plot(x_points_4nodes, T_points_4nodes, 'LineWidth',2); grid on; xlim([x_points_4nodes(1) x_points_4nodes(end)]);
xlabel('$x \left[ mm \right]$', 'FontSize',18); ylabel('Temperature $ \left[ K \right]$', 'FontSize',18);
xticks([0 distance_1, distance_2, distance_3, distance_4, distance_5, distance_6, distance_7]);
xticklabels({'0','$l_1$','$\frac{layer_2}{2}$', '$l_2$', '$l_3$','$\frac{layer_4}{2}$', '$l_4$', '$l_5$'});
ax = gca; ax.FontSize = 15;



%% FUNCTIONS EX1 

function [accelerations] = state_acceleration(state_fun, state_values)
% Function that computes the accelerations of the states at each time step
% of the solution given the state function x_dot =
% state_fun(t,state_value);
% In this case we already assume that accelerations in vector x_dot are in
% positions 2 and 4
accelerations = zeros(length(state_values), 2);

for j = 1:length(state_values)
    state = state_values(j,:);
    x_dot = state_fun(state);
    accelerations(j,1) = x_dot(2);
    accelerations(j,2) = x_dot(4);
end


end



function [minimum, row, column] = matrix_minimum(A)
% function that computes and locate the minimum value inside a matrix
[min_column, row_index] = min(A);
[minimum, column] = min(min_column);
row = row_index(column);

end

%% FUNCTIONS EX2

function [x_dot, useful_parameters] = hydraulic_network(t, x, parameters)

% function that expresses the state equation for the hydraulic network.
% The state is assumed as [V_accumulator, x_actuator, v_actuator, V_tank].
% Parameters are to be provided as [rho, accumulator, delivery, distributor, actuator, rline, tank]

V_acc = x(1);
x_act = x(2);
v_act = x(3);

rho = parameters{1};
accumulator = parameters{2};
delivery = parameters{3};
distributor = parameters{4};
actuator = parameters{5};
rline = parameters{6};
tank = parameters{7};

% check on the x value
if x_act < 0
    x_act = 0;
end
if x_act <= 0 && v_act < 0
    x_act = 0;
    v_act = 0;
end
if x_act > actuator.x_max 
    x_act = actuator.x_max;
end
if x_act >= actuator.x_max && v_act > 0
    x_act = actuator.x_max;
    v_act = 0;
end

% define the command z:
if t < 1
    z = 0;
elseif t > 1 && t < 1.5
    z = 2*(t-1);
elseif t> 1.5
    z = 1;
end

% area of the distributor depending on z
alpha = 2*acos(1 - 2*z);  % angle of the opening
A_distributor = distributor.r^2/2*( alpha - sin(alpha) ); % area of the distributor

% Before writing the state equation we need to properly set up all flow
% rates and pressures.
% It will be considered that z is always >= 0

A_actuator4 = pi*actuator.R_cyl^2;
A_actuator5 = pi*actuator.R_cyl^2 - pi*actuator.R_rod^2;  % area of the 2 actuator chambers

Q4 = sign(z)*A_actuator4*v_act;   % flow rate depends on actuator speed.
Q5 = sign(z)*A_actuator5*v_act;   % They are both set to zero when z == 0.
Q3 = Q4;
Q6 = Q5;
Q2 = Q3;
Q1 = Q2;
Q_acc = Q1;
Q7 = Q6;
Q_tank = Q7;

A_delivery = delivery.A;   % delivery pipe area.
A_return   = rline.A;      % return pipe area.

% pressures, connected by head and minor losses. Before the distributor 
% computed started from the accumulator pressure. After the distributor
% we start from the tank pressure, assumed constant.

P_acc = accumulator.p0*( accumulator.V0 / ( accumulator.V0 + V_acc ) )^accumulator.gamma; % we have a plus on the denominator because as the flow exits the pressure decreases.
P1 = P_acc  - 0.5*delivery.k_acc*rho*Q1*abs(Q1)/(A_delivery)^2; % minor loss, accumulator exit.
P2 = P1     - 0.5*delivery.k_cv*rho*Q2*abs(Q2)/(A_delivery)^2;  % minor loss, check-valve.
P3 = P2     - 0.5*delivery.f*rho*delivery.L/delivery.D*Q3*abs(Q3)/(A_delivery)^2; % head loss: 2-3 pipe.
P7 = tank.p + 0.5*rline.kt*rho*Q7*abs(Q7)/(A_return)^2;                           % minor loss, tank inlet.
P6 = P7     + 0.5*rline.f*rho*rline.L/rline.D*Q6*abs(Q6)/(A_return)^2;            % head loss, pipe 6-7.


if z == 0 || A_distributor < 1e-7   % check on the distributor for P4: if closed P4 grants equilibrium with P5, when open depends on incoming flow.
    P5 = tank.p;                    % need also to set a proper tolerance on A_distributor to avoid numerical issues; If omitted delivers weird results.
    P4 = ( P5*A_actuator5 + actuator.load(x_act) )/A_actuator4;
elseif z > 0       
    P4 = P3     - 0.5*distributor.k*rho*Q4*abs(Q4)/(A_distributor)^2;                  % minor loss, distributor on delivery line.
    P5 = P6     + 0.5*distributor.k*rho*Q5*abs(Q5)/(A_distributor)^2;                  % minor loss, distributor on return line.
end

% now we can define the equation for the state variables as x_dot = f(t,x)
x_dot = [ Q_acc;...
          v_act;...
          ( P4*A_actuator4 - P5*A_actuator5 - actuator.load(x_act) )/actuator.m ;...
          Q_tank];
if x_act >= actuator.x_max % needed to control the actuator velocity after the stroke occurs
    x_dot(2) = 0;
    x_dot(3) = 0;
end

useful_parameters = [ Q_acc, Q1, Q2, Q3, Q4, Q5, Q6, Q7, P_acc, P1, P2, P3, P4, P5, P6, P7 ];


end


function [value, isterminal, direction] = actuator_stroke(~, x, x_max)
value = x_max - x(2);
isterminal = 1;
direction = 0;
end

%% FUNCTIONS EX4

function [dT] = simulation_2nodes(t,x,parameters)
% simulates behaviour in a 2 node model. x = [T1,T2], parameters =
% {density,conductivity, specific_h,layer_l, T0, Tf}.

density = parameters{1};
conductivity = parameters{2};
specific_h = parameters{3};
layer_l = parameters{4};
T0 = parameters{5};
Tf = parameters{6};

C1 = density(2)*specific_h(2)*layer_l(2);   % capacity layer 2
C2 = density(4)*specific_h(4)*layer_l(4);   % capacity layer 4

if t < 1
    T_i = T0 + (Tf - T0)*t;
else
    T_i = Tf;
end

% thermal resistances
R1 = layer_l(1)/conductivity(1);   % resistance layer 1
R2 = layer_l(2)/conductivity(2);   % resistance layer 2
R3 = layer_l(3)/conductivity(3);   % resistance layer 3
R4 = layer_l(4)/conductivity(4);   % resistance layer 4
R5 = layer_l(5)/conductivity(5);   % resistance layer 5

T1 = x(1);  % temperature node at layer 2
T2 = x(2);  % temperature node at layer 4

dT = [ ( (R1 + R2/2)^(-1)*(T_i - T1) - (R2/2 + R3 * R4/2)^(-1)*(T1 - T2) )/C1;
       ( (R2/2 + R3 + R4/2)^(-1)*(T1 - T2) - (R4/2 + R5)^(-1)*(T2 - T0) )/C2 ];
end



function [dT] = simulation_4nodes(t,x,parameters)
% simulates behaviour in a 4 node model. x = [T1,T2], parameters =
% {density,conductivity, specific_h,layer_l, T0, Tf}.

density = parameters{1};
conductivity = parameters{2};
specific_h = parameters{3};
layer_l = parameters{4};
T0 = parameters{5};
Tf = parameters{6};

C1 = density(2)*specific_h(2)*layer_l(2);   % capacity layer 2
C2 = density(4)*specific_h(4)*layer_l(4);   % capacity layer 4

if t < 1
    T_i = T0 + (Tf - T0)*t;
else
    T_i = Tf;
end

% thermal resistances
R1 = layer_l(1)/conductivity(1);   % resistance layer 1
R2 = layer_l(2)/conductivity(2);   % resistance layer 2
R3 = layer_l(3)/conductivity(3);   % resistance layer 3
R4 = layer_l(4)/conductivity(4);   % resistance layer 4
R5 = layer_l(5)/conductivity(5);   % resistance layer 5

T1 = x(1);  % temperature 1-st node at layer 2
T2 = x(2);  % temperature 2-nd node at layer 2
T3 = x(3);  % temperature 1-st node at layer 4
T4 = x(4);  % temperature 2-nd node at layer 4


dT = [ ( (R1 + R2/3)^(-1)*(T_i - T1) - (R2/3)^(-1)*(T1 - T2) )/(C1/2);
       ( (R2/3)^(-1)*(T1 - T2) - (R2/3 + R3 + R4/3)^(-1)*(T2 - T3) )/(C1/2);
       ( (R2/3 + R3 + R4/3)^(-1)*(T2 - T3) - (R4/3)^(-1)*(T3 - T4) )/(C2/2);
       ( (R4/3)^(-1)*(T3 - T4) - (R4/3 + R5)^(-1)*(T4 - T0) )/(C2/2)  ];
end


function [T_i] = imposed_temperature(t_span, T0, Tf)
T_i = zeros(length(t_span),1);
for i = 1:length(t_span)
    if t_span(i) < 1
        T_i(i) = T0 + (Tf - T0)*t_span(i);
    else
        T_i(i) = Tf;
    end
end

end