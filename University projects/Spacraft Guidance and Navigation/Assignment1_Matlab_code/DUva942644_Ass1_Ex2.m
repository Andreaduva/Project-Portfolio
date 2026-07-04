%% ASSIGNMENT 1 EX 2

% initialize kernels and other settings
close all;
clear; clc;
cspice_kclear();

% change plots appearance 
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');
% set parallel computing pool
pp = gcp;
if pp.Connected ~= 1
   parpool;
end

% set constants used during the exercise
mass_ratio = 1.21506683e-2; % Earth-Moon mass parameter            [-]
mass_sun = 3.28900541e5;    % scaled mass of the Sun               [-]
rho = 3.88811143e2;         % scaled Sun-(Earth+Moon) distance     [-]
omega_sun = -9.25195985e-1; % scaled angular velocity of the Sun   [-]
l_EM = 3.84405000e08;       % Earth-Moon distance                  [m]
omega_EM = 2.66186135e-6;   % Earth-Moon angular velocity          [s^-1]
Re = 6378;                  % Mean Earth's radius                  [km]
Rm = 1738;                  % Mean Moon's radius                   [km]
hi = 167;                   % Altitude of departure orbit          [km]
hf = 100;                   % Altitude of arrival orbit            [km]
DU = 3.8440500e08;          % Distance unit                        [m]
TU = 4.34811305;            % Time unit                            [days]
VU = 1.02323281e3;          % Speed unit                           [m/s]

% retrieve other useful parameters
r_i = (Re + hi)*1e3/DU;           % scaled height of the Earth's parking orbit [-]
v_i = sqrt((1 - mass_ratio)/r_i); % scaled velocity of the Earth's orbit       [-]
r_f = (Rm + hf)*1e3/DU;           % scaled height of the Moon's parking orbit  [-]
v_f = sqrt(mass_ratio/r_f);       % scaled velocity of the Moon's orbit        [-]


%% EX 2 POINT 1

% set parameters for the initial conditions:
alpha = 1.5*pi;
beta = 1.41;
ti = 0;
delta = 7;
parameters = [mass_ratio; mass_sun; rho; omega_sun];

% compute them in the form useful for the integration:
[x0,t_i,t_f] = initial_condition(alpha, beta, ti, delta, r_i, mass_ratio);

% propagate the obtained initial conditions
opt_PBRFBP = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
[tt,xx_rot] = ode113(@(t,x) PBRFBP_ode(t,x,parameters), [t_i t_f],x0,opt_PBRFBP);

% plot the orbit in the rotating Earth-Moon frame:
figure();
grid on;
hold on;
plot(xx_rot(:,1),xx_rot(:,2),'LineWidth',2);
plot(-mass_ratio,0,'ok','MarkerFaceColor','k','MarkerSize',7); % Earth
plot(1 - mass_ratio,0,'ob','MarkerFaceColor','b','MarkerSize',7); % Moon
% Sun position at t_i and t_f:
Sun_pos_fun = @(t) rho*[ cos(omega_sun.*t);
                     sin(omega_sun.*t)];
Sun_pos = Sun_pos_fun([t_i t_f]);
quiver(0,0,Sun_pos(1,1)/rho,Sun_pos(1,2)/rho,'Color','#EDB120','LineWidth',3);
quiver(0,0,3*Sun_pos(2,1)/rho,3*Sun_pos(2,2)/rho,'Color',"#D95319",'LineWidth',3);
xlabel('$x [-]$'); ylabel('$y [-]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('Trajectory','$m_1$','$m_2$','Sun at $t_i$','Sun at $t_f$','FontSize',15,'Location','northwest');

% plot the orbit in the earth centered inertial frame:
xx_in = rotating2inertial2D(tt,xx_rot,mass_ratio);
Moon_rot = [1-mass_ratio,0,0,0];
Moon_state_tf = rotating2inertial2D(t_i,Moon_rot,mass_ratio);
figure();
grid on;
hold on;
plot(xx_in(:,1),xx_in(:,2),'LineWidth',2);
plot(0,0,'ob','MarkerFaceColor','b','MarkerSize',7); % Earth state
plot(cos(tt(:)),sin(tt(:)),'LineWidth',2);  % Moon Orbit
xlabel('$X [-]$'); ylabel('$Y [-]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('Trajectory','Earth','Moon orbit','FontSize',15,'Location','northwest');
% comparison between numerical STM and variational one:
[x_f_var, STM_var] = PBRFBP_STM_variational(t_i, t_f, x0, parameters);
%% EX 2 POINT 2: SIMPLE SHOOTING

% 2.2a: run optimizer without providing derivatives
parameters_ss = [parameters; r_i; r_f];        % parameters for the Simple Shooting function

% define upper and lower bounds for the variables
x_lb = - r_i - mass_ratio;
x_ub =   r_i - mass_ratio;
y_lb = - r_i;
y_ub =   r_i;
t_i_ub = 2*pi/abs(omega_sun);
v0_ub = sqrt( 2*(1 - mass_ratio)/r_i );
Vx_lb = -(v0_ub - r_i);
Vx_ub =  (v0_ub - r_i);
Vy_lb = Vx_lb;
Vy_ub = Vx_ub;
lb_vect = [x_lb;y_lb;Vx_lb;Vy_lb;0;0];
ub_vect = [x_ub;y_ub;Vx_ub;Vy_ub;t_i_ub;23*TU];
% linear inequality constraints: t_f > t_i
A_ineq = [0 0 0  0  1 -1];
b_ineq = 0;
% nonlinear constraints function
nonlincon_ss = @(y) Earth_Moon_nonlincon_ss(y,parameters_ss);
% cost function
cost_fun_ss = @(y) Earth_Moon_cost_ss(y,parameters_ss);
% first guess:
y0_ss = [x0;t_i;t_f];
% call fmincon
opt_ss_num = optimoptions('fmincon','Display','iter-detailed','MaxFunctionEvaluations',7000,'Algorithm','active-set','ConstraintTolerance',1e-10, ...
    'StepTolerance',2e-14,'FunctionTolerance',1e-14,'UseParallel',true);
%'FiniteDifferenceType','central'
tic;
[y_sol_a,fval_a,flag_a,~] = fmincon(cost_fun_ss,y0_ss,A_ineq,b_ineq,[],[],lb_vect,ub_vect,nonlincon_ss,opt_ss_num);
toc;

% propagate the obtained results:
opt_PBRFBP = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
[~,xx_opt_a] = ode113(@(t,x) PBRFBP_ode(t,x,parameters), [y_sol_a(5) y_sol_a(6)],y_sol_a(1:4),opt_PBRFBP);
% plot the orbit in the rotating Earth-Moon frame:
figure();
grid on;
hold on;
plot(xx_opt_a(:,1),xx_opt_a(:,2),'LineWidth',2);
plot(-mass_ratio,0,'ok','MarkerFaceColor','k','MarkerSize',7); % Earth
plot(1 - mass_ratio,0,'ob','MarkerFaceColor','b','MarkerSize',7); % Moon
Sun_pos_a = Sun_pos_fun([y_sol_a(5) y_sol_a(6)]);
quiver(0,0,Sun_pos_a(1,1)/rho,Sun_pos_a(1,2)/rho,'Color','#EDB120','LineWidth',3);
quiver(0,0,Sun_pos_a(2,1)/rho,Sun_pos_a(2,2)/rho,'Color',"#D95319",'LineWidth',3);
xlabel('$x [-]$'); ylabel('$y [-]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('Trajectory','$m_1$','$m_2$','Sun at $t_i$','Sun at $t_f$','FontSize',15,'Location','northwest');



% 2.2b : provide gradients to the solver

opt_ss_var = optimoptions('fmincon','Display','iter-detailed','MaxFunctionEvaluations',7000,'Algorithm','active-set','ConstraintTolerance',1e-10, ...
    'StepTolerance',2e-14,'FunctionTolerance',1e-14,'SpecifyConstraintGradient',true,'SpecifyObjectiveGradient',true);
% tic;
% [y_sol_b,fval_b,flag_b,~] = fmincon(cost_fun,[y0(1:4);y_sol_a(5);y0(6)],A_ineq,b_ineq,[],[],lb_vect,ub_vect,nonlincon,opt);
% toc; 

tic;
[y_sol_b,fval_b,flag_b,~] = fmincon(cost_fun_ss,y0_ss,A_ineq,b_ineq,[],[],lb_vect,ub_vect,nonlincon_ss,opt_ss_var);
toc;

% propagate the obtained results:
opt_PBRFBP = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
[~,xx_opt_b] = ode113(@(t,x) PBRFBP_ode(t,x,parameters), [y_sol_b(5) y_sol_b(6)],y_sol_b(1:4),opt_PBRFBP);
% plot the orbit in the rotating Earth-Moon frame:
figure();
grid on;
hold on;
plot(xx_opt_b(:,1),xx_opt_b(:,2),'LineWidth',2);
plot(-mass_ratio,0,'ok','MarkerFaceColor','k','MarkerSize',7); % Earth
plot(1 - mass_ratio,0,'ob','MarkerFaceColor','b','MarkerSize',7); % Moon
Sun_pos_b = Sun_pos_fun([y_sol_a(5) y_sol_a(6)]);
quiver(0,0,Sun_pos_b(1,1)/rho,Sun_pos_b(1,2)/rho,'Color','#EDB120','LineWidth',3);
quiver(0,0,Sun_pos_b(2,1)/rho,Sun_pos_b(2,2)/rho,'Color',"#D95319",'LineWidth',3);
xlabel('$x [-]$'); ylabel('$y [-]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('Trajectory','$m_1$','$m_2$','Sun at $t_i$','Sun at $t_f$','FontSize',15,'Location','northwest');

%% POINT 3: MULTIPLE SHOOTING

% first generate initial guesses for vector y using the starting point x0 in Point 1
N = 4;
t_vec = linspace(t_i,t_f,N);      % time grid
x0_MS = zeros(4*N,1);
x0_MS(1:4) = x0;
% compute flow at each time point
for i = 2:N
    flow_im1 = PBRFBP_flow(t_vec(i-1),t_vec(i),x0_MS(4*(i-1)-3:4*(i-1)),parameters);  % state at node i from propagation of node i-1
    x0_MS(4*i-3:4*i) = flow_im1;
end

% plot initial guess
figure();
grid on;
hold on;
plot(1 - mass_ratio,0,'ob','MarkerFaceColor','b','MarkerSize',7); % Moon
for i = 1:N-1
    opt_PBRFBP = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
    [~,xx_leg] = ode113(@(t,x) PBRFBP_ode(t,x,parameters), [t_vec(i) t_vec(i+1)],x0_MS(4*i-3:4*i),opt_PBRFBP);
    plot(xx_leg(:,1),xx_leg(:,2),'Color',"#0072BD",'LineWidth',1.5);
    plot(x0_MS(4*i-3),x0_MS(4*i-2),'or','MarkerFaceColor','r','MarkerSize',6)
end
plot(x0_MS(end-3),x0_MS(end-2),'or','MarkerFaceColor','r','MarkerSize',6)
xlabel('$x [-]$'); ylabel('$y [-]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('$m_2$','Trajectory','MS nodes','FontSize',15,'Location','northwest');

% set up quantities for minimization process
parameters_ms = [parameters_ss;Re/DU;Rm/DU];

% nonlinear constraints function
nonlincon_ms = @(y) Earth_Moon_nonlincon_ms(y,parameters_ms);
% cost function
cost_fun_ms = @(y) Earth_Moon_cost_ms(y,parameters_ms);
% first guess:
y0_ms = [x0_MS;t_i;t_f];
% call fmincon
opt_ms = optimoptions('fmincon','Display','iter-detailed','MaxFunctionEvaluations',7000,'Algorithm','active-set','ConstraintTolerance',1e-10, ...
    'StepTolerance',2e-14,'FunctionTolerance',1e-14,'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true);
tic;
[y_sol_ms,fval_ms,flag_ms,~] = fmincon(cost_fun_ms,y0_ms,[],[],[],[],[],[],nonlincon_ms,opt_ms);
toc;
t_i_sol_ms = y_sol_ms(end-1);
t_f_sol_ms = y_sol_ms(end);

% plot solution of multiple shooting
t_vec_sol = linspace(t_i_sol_ms,t_f_sol_ms,N);      % time grid
figure();
grid on;
hold on;
plot(1 - mass_ratio,0,'ob','MarkerFaceColor','b','MarkerSize',7); % Moon
for i = 1:N-1
    opt_PBRFBP = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
    [~,xx_leg] = ode113(@(t,x) PBRFBP_ode(t,x,parameters), [t_vec_sol(i) t_vec_sol(i+1)],y_sol_ms(4*i-3:4*i),opt_PBRFBP);
    plot(xx_leg(:,1),xx_leg(:,2),'Color',"#0072BD",'LineWidth',1.5);
    plot(y_sol_ms(4*i-3),y_sol_ms(4*i-2),'or','MarkerFaceColor','r','MarkerSize',6)
end
plot(y_sol_ms(end-5),y_sol_ms(end-4),'or','MarkerFaceColor','r','MarkerSize',6)
xlabel('$x [-]$'); ylabel('$y [-]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('$m_2$','Trajectory','MS nodes','FontSize',15,'Location','northeast');

% check quantities
flow_ss1 = PBRFBP_flow(y_sol_a(5),y_sol_a(6),y_sol_a(1:4),parameters);
flow_ss2 = PBRFBP_flow(y_sol_b(5),y_sol_b(6),y_sol_b(1:4),parameters);
xN_sol = y_sol_ms(end-5:end-2);
flow_first_guess = PBRFBP_flow(t_i,t_f,x0,parameters);
xN_first_guess = y0_ms(end-5:end-2);
DeltaV_first_guess = Earth_Moon_cost_ss(y0_ss,parameters_ss);
DeltaV_first_ms = Earth_Moon_cost_ms(y0_ms,parameters_ms);

% comparing final states:
figure();
grid on;
hold on;
h1=plot(1 - mass_ratio,0,'ob','MarkerFaceColor','b','MarkerSize',7); % Moon
h2=viscircles([1-mass_ratio 0],r_f,'Color',[0.5 0.5 0.5],'LineWidth',2); % Moon parking orbit
h3=scatter(y_sol_ms(end-5),y_sol_ms(end-4),55,'o','filled');
h4=scatter(xx_opt_b(end,1),xx_opt_b(end,2),55,'o','filled');
xlabel('$x [-]$'); ylabel('$y [-]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend([h1,h2,h3,h4],{'Moon','Parking orbit','MS result','SS result'},'FontSize',15,'Location','northeast');
%% POINT 1 FUNCTIONS

function [x0,ti,tf] = initial_condition(alpha, beta, ti, delta, r0, mu)
% Function that given four scalar parameters can reconstruct the initial state in the classical
% expression. All quantities are nondimensional.
% INPUT:
% alpha = angle defined on Earth circular parking orbit [-]
% beta  = initial-to-circular velocity ratio            [-]
% ti    = initial time                                  [-]
% delta = transfer duration                             [-]
% r0    = orbit's height at initial time                [-]
% mu    = mass ratio of the system                      [-]
%
% OUTPUT
% x0 = 4x1 vector of initial conditions (r0;v0)         [-]
% ti = initial time                                     [-]
% tf = final time of the transfer                       [-]

% initialize parameters
x0 = zeros(4,1);
v0 = beta*sqrt( (1 - mu)/r0 );

x0(1) = r0*cospi(alpha/pi) - mu;     % x0
x0(2) = r0*sinpi(alpha/pi);          % y0
x0(3) = -(v0 - r0)*sinpi(alpha/pi);  % x_dot0
x0(4) = (v0 - r0)*cospi(alpha/pi);   % y_dot0
tf    = ti + delta;
end


function [x_in] = rotating2inertial2D(t,x_rot,mu)
% Function that converts position and velocity coordinates from the rotating frame to an inertial
% Earth centered one.
%
% INPUT:
% t     = vector of times corresponding to the given coordinates.   [-]
% x_rot = matrix containing the set of coordinates in rotating frame; each row corresponds to a specific epoch t.
%         In this case the columns are (x,y,Vx,Vy).                 [-]
% mu    = mass parameter of the inner P1-P2 system.                 [-]
% OUTPUT:
% x_in  = matrix containing the set of coordinates in the Earth centered inertial frame.  [-]

x_in = zeros(length(t),4);
x_in(:,1) = (x_rot(:,1) + mu).*cos(t(:)) - x_rot(:,2).*sin(t(:));  % x coordinate
x_in(:,2) = (x_rot(:,1) + mu).*sin(t(:)) + x_rot(:,2).*cos(t(:));  % y coordinate
x_in(:,3) = (x_rot(:,3) - x_rot(:,2)).*cos(t(:)) - (x_rot(:,4) + x_rot(:,1) + mu).*sin(t(:));  % Vx coordinate
x_in(:,4) = (x_rot(:,3) - x_rot(:,2)).*sin(t(:)) + (x_rot(:,4) + x_rot(:,1) + mu).*cos(t(:));  % Vy coordinate
end



function [x_dot] = PBRFBP_ode(t,state,parameters)
% Functions that expresses the ode governing the motions of the S/C in the normalized PBRFBP, written in canonical form x_dot = f(x,t)
% All quantities are nondimensional.
% INPUT:
% t = current time                                [-]
% state = 4x1 state vector (x,y,Vx,Vy)            [-]
% parameters = vector containing numerical values of the following:
%         mu = scalar governing the P1-P2 motion          [-]
%         mass_S = scaled mass of P3 (Sun)                [-]
%         rho = scaled P3-(P1-P2) distance                [-]
%         omega_S = scaled angular velocity of P3 (Sun)   [-]
%
% OUTPUT:
% x_dot = 4x1 vector of state derivatives (Vx,Vy,acc_x,acc_y)  [-]

% check size of the parameters vector:
if length(parameters) ~= 4
    error('\nThe parameters vector passed to the ode function must have 4 elements\n');
end

mu = parameters(1);         % Earth-Moon mass parameter            [-]
mass_3b = parameters(2);    % scaled mass of the Sun (third body)  [-]
rho = parameters(3);        % scaled Sun-(Earth+Moon) distance     [-]
omega_3b = parameters(4);   % scaled angular velocity of the Sun   [-]


% set useful quantities
x_dot = zeros(4,1);
x =  state(1);                           % S/C x position
y =  state(2);                           % S/C y position
Vx = state(3);                           % S/C x velocity
Vy = state(4);                           % S/C y velocity
xS = rho*cos(omega_3b*t);                % Sun x position
yS = rho*sin(omega_3b*t);                % Sun y position
r1 = sqrt( (x + mu).^2 + y.^2);          % S/C-Earth distance
r2 = sqrt( (x + mu - 1).^2 + y.^2);      % S/C-Moon distance
r3 = sqrt( (x - xS).^2 + (y - yS).^2);   % S/C-Sun distance

% write state_derivative = f(x,t)
x_dot(1) = Vx;
x_dot(2) = Vy;
x_dot(3) = 2*Vy + x - (1 - mu)*(mu + x)./(r1.^3) + mu*(1 - mu - x)./(r2.^3) - mass_3b*(x - xS)./(r3.^3) - mass_3b/(rho^2)*cos(omega_3b*t);
x_dot(4) = -2*Vx + y - y*(1 - mu)./(r1.^3) - mu*y./(r2.^3) - mass_3b*(y - yS)./(r3.^3) - mass_3b/(rho^2)*sin(omega_3b*t);
end



function [state_f]  = PBRFBP_flow(t_i,t_f,state0,parameters)
% Functions that expresses the flow of the ode governing the motions of the S/C in the normalized PBRFBP, written in canonical form x_dot = f(x,t)
% All quantities are nondimensional.
% INPUT:
% t_i         = initial time of integration      [-]
% t_f         = final time of integration        [-]
% state0 = 4x1 initial state vector (x0,y0,Vx0,Vy0)       [-]
% parameters = vector containing numerical values of the following:
%         mu = mass ratio governing the P1-P2 motion      [-]
%         mass_S = scaled mass of P3 (Sun)                [-]
%         rho = scaled P3-(P1-P2) distance                [-]
%         omega_S = scaled angular velocity of P3 (Sun)   [-]
%
% OUTPUT:
% state_f = 4x1 vector of final state (x_f,y_f,Vx_f,Vy_f) [-]


% check size of the parameters vector:
if length(parameters) ~= 4
    error('\nThe parameters vector passed to the flow function must have 4 elements\n');
end

% propagate initial S/C condition to the final time:
opt = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
[~,xx] = ode113(@(t,x) PBRFBP_ode(t,x,parameters), [t_i t_f],state0,opt);
% retrieve the final state
state_f = xx(end,:)';
end


function [dXdT] = PBRFBP_ode_variational(t,state,parameters)
% Function that expresses the ode for the PBRCTBP and its Jacobiam matrix when requested in order
% to obtain the flow at a given time and the corresponding STM
% INPUT:
% t          = current time values                                           [-]
% state      = 20x1 vector containing the current state (r;v;STM components) [-]
% parameters = vector of parameters built as: 
%         mu = mass ratio governing the P1-P2 motion       [-]
%         mass_S = scaled mass of P3 (Sun)                 [-]
%         rho = scaled P3-(P1-P2) distance                 [-]
%         omega_S = scaled angular velocity of P3 (Sun)    [-]
%
% OUTPUT:
% dXdT  = 20x1 vector for the extended state derivative (v;v_dot;STM_dot)  [-]

% check size of the parameters vector:
if length(parameters) ~= 4
    error('\nThe parameters vector passed to the variational equation must have 4 elements\n');
end
% retrieve parameters
PHI = reshape(state(5:end), [4 4]);
mu = parameters(1);         % Earth-Moon mass parameter            [-]
mass_3b = parameters(2);    % scaled mass of the Sun (third body)  [-]
rho = parameters(3);        % scaled Sun-(Earth+Moon) distance     [-]
omega_3b = parameters(4);   % scaled angular velocity of the Sun   [-]

% set useful quantities
DxDt = zeros(4,1);
x =  state(1);                           % S/C x position
y =  state(2);                           % S/C y position
Vx = state(3);                           % S/C x velocity
Vy = state(4);                           % S/C y velocity
xS = rho*cos(omega_3b*t);                % Sun x position
yS = rho*sin(omega_3b*t);                % Sun y position
r1 = sqrt( (x + mu).^2 + y.^2);          % S/C-Earth distance
r2 = sqrt( (x + mu - 1).^2 + y.^2);      % S/C-Moon distance
r3 = sqrt( (x - xS).^2 + (y - yS).^2);   % S/C-Sun distance

% write state_derivative = f(x,t)
DxDt(1) = Vx;
DxDt(2) = Vy;
DxDt(3) = 2*Vy + x - (1 - mu)*(mu + x)./(r1.^3) + mu*(1 - mu - x)./(r2.^3) - mass_3b*(x - xS)./(r3.^3) - mass_3b/(rho^2)*cos(omega_3b*t);
DxDt(4) = -2*Vx + y - y*(1 - mu)./(r1.^3) - mu*y./(r2.^3) - mass_3b*(y - yS)./(r3.^3) - mass_3b/(rho^2)*sin(omega_3b*t);

% we compute also the Jacobian of the dynamics equation. The non trivial terms
% are:
OMxx = 1 - (1 - mu)/(r1.^3) + 3*(1 - mu)*((x + mu).^2)./(r1.^5) - mu/(r2.^3) + 3*mu*((x + mu - 1).^2)./(r2.^5)-mass_3b/(r3.^3) + 3*mass_3b*(x - xS).^2./(r3.^5);
OMyy = 1 - (1 - mu)/(r1.^3) + 3*(1 - mu)*(y.^2)./(r1.^5) - mu/(r2.^3) + 3*mu*(y.^2)./(r2.^5) - mass_3b./(r3.^3) + 3*mass_3b.*(y - yS).^2./(r3.^5);
OMxy = 3*(1 - mu)*(x + mu).*y./(r1.^5) + 3*mu*(x + mu - 1).*y./(r2.^5) + 3*mass_3b*(x - xS).*(y - yS)./(r3.^5);   % equal to OMyx

% the jacobian matrix can now be assembled
A = [zeros(2),eye(2);
     OMxx, OMxy,  0, 2;
     OMxy, OMyy, -2, 0];  % 4x4 matrix

% derivative of the STM:
PHI_dot = A*PHI;
% get back to a vectorial expression of the augmented state
dXdT = [DxDt; PHI_dot(:)];
end


function [x_f, STM] = PBRFBP_STM_variational(t_i,t_f,x0,parameters)
% Function that computes the flow of the ODE and the STM through finite differences
% INPUT:
% t_i         = initial time of integration    [-]
% t_f         = final time of integration      [-]
% x0          = 4x1 vector containing the initial state state (r0;v0)   [-]
% parameters  = vector of parameters built as: 
%          mu = mass ratio governing the P1-P2 motion       [-]
%          mass_S = scaled mass of P3 (Sun)                 [-]
%          rho = scaled P3-(P1-P2) distance                 [-]
%          omega_S = scaled angular velocity of P3 (Sun)    [-]
%
% OUTPUT:
% x_f  = flow of the ode at the final time. 4x1 vector (r_f; v_f)   [-]
% STM  = State Transition Matrix at the final time, 4x4 matrix      [-]
%
% External functions called: PBRFBP_ode_variational


% initialize initial conditions:
state0 = [x0; reshape( eye(4), [16 1])];
% numerical integration:
opt = odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
[~,xx] = ode113(@(t,x) PBRFBP_ode_variational(t,x,parameters), [t_i, t_f], state0, opt);
x_f = xx(end,1:4)'; 
STM = reshape(xx(end,5:end), [4 4]);
end

%% POINT 2 FUNCTIONS

function[F,G] = Earth_Moon_cost_ss(state1,parameters)
% Function that computes the total cost associated to the two impulse Earth-Moon transfer specified
% by the values contained in the vector y. The initial and final delta velocity are computed with
% respect to prescribed and fixed circular orbits. Single shooting method considered.
%
% INPUT:
% state1 = 6x1 augmented state vector built as: (x0, y0, Vx0, Vy0, t_i, t_f). [-]
% parameters = vector of parameters built as: 
%         mu = mass ratio governing the P1-P2 motion       [-]
%         mass_S = scaled mass of P3 (Sun)                 [-]
%         rho = scaled P3-(P1-P2) distance                 [-]
%         omega_S = scaled angular velocity of P3 (Sun)    [-]
%         r_E = height of Earth's parking orbit            [-]
%         r_M = height of Moon's parking orbit             [-]
% OUTPUT:
% F    = scalar cost index associated to the total cost of the transfer defined as
%          |v_i - v0| + |v_f - velocity_flow(t_i;t_f)|                               [-]
% G    = Gradient of the cost funtion with respect to the NLP variables. 6x1 vector. [-]
%        Computed only when required.

% check on the size of the parameters vector 
if length(parameters) ~= 6
    error('\n The parameters vector passed to the objective function\nin the simple shooting case must have 6 elements\n');
end

% define variables for better reading
x1 = state1(1);
y1 = state1(2);
Vx1 = state1(3);
Vy1 = state1(4);
t_i = state1(5);
t_f = state1(6);
mu = parameters(1);        % mass ratio
r_E = parameters(5);       % height of Earth's parking orbit
r_M = parameters(6);       % height of Moon's target orbit

% compute the state at final time:
if nargout == 1
    state2 = PBRFBP_flow(t_i,t_f, state1(1:4),parameters(1:4));  % only state at final time needed
else
    [state2,STM] = PBRFBP_STM_variational(t_i,t_f,state1(1:4),parameters(1:4)); % to compute the gradient G we need also the STM
end
x_2  = state2(1);
y_2  = state2(2);
Vx_2 = state2(3);
Vy_2 = state2(4);
% compute the initial and final DeltaVs
DeltaV_i = sqrt( (Vx1 - y1).^2 + (Vy1 + x1 + mu).^2 ) - sqrt( (1 - mu)/r_E );
DeltaV_f = sqrt( (Vx_2 - y_2).^2 + (Vy_2 + x_2 + mu - 1).^2 ) - sqrt(mu/r_M);
% overall cost is then:
F = DeltaV_i + DeltaV_f;


% computing the gradient:
if nargout > 1
    G = zeros(6,1);     % initialize gradient
    Dv1_Dx1 = 1/sqrt( (Vx1 - y1).^2 + (Vy1 + x1 + mu).^2 ).*[ Vy1 + x1 + mu;
                                                              y1 - Vx1;
                                                              Vx1 - y1;
                                                              Vy1 + x1 + mu];
    Dv2_Dx2 = 1/sqrt( (Vx_2 - y_2).^2 + (Vy_2 + x_2 + mu - 1).^2 ).*[ Vy_2 + x_2 + mu - 1;
                                                                      Vy_2 - Vx_2;
                                                                      Vx_2 - y_2;
                                                                      Vy_2 + x_2 + mu - 1];
    G(1:4) = Dv1_Dx1 + STM'*Dv2_Dx2;
    G(5) = - PBRFBP_ode(t_i,state1(1:4),parameters(1:4))'*STM'*Dv2_Dx2;
    G(6) = + PBRFBP_ode(t_f,state2,parameters(1:4))'*Dv2_Dx2;
end
end


function [c,c_eq,Gc,Gc_eq] = Earth_Moon_nonlincon_ss(state1,parameters)
% Function expressing the nonlinear constraints for the Earth-Moon transfer optimization using the
% simple shooting technique. The 4 scalar constraints to be enforced concerns the physical distance
% and the velocity direction at both initial and final time. Single shooting method considered.
%
% INPUT:
% state1 = 6x1 augmented state vector built as: (x0, y0, Vx0, Vy0, t_i, t_f). [-]
% parameters = vector of parameters built as: 
%         mu = mass ratio governing the P1-P2 motion       [-]
%         mass_S = scaled mass of P3 (Sun)                 [-]
%         rho = scaled P3-(P1-P2) distance                 [-]
%         omega_S = scaled angular velocity of P3 (Sun)    [-]
%         r_E = height of Earth's parking orbit            [-]
%         r_M = height of Moon's parking orbit             [-]
% OUTPUT:
% c     = nonlinear inequality constraints. Not defined in this problem.   [-] 
% c_eq  = 4x1 vectorial funtion expressing nonlinear equality constraints. [-]
% Gc    = gradient of nonlinear inequality constraints, still not defined  [-]
% Gc_eq = gradient of equality constraints, provided only if requested.    [-]


% check on the size of the parameters vector 
if length(parameters) ~= 6
    error('\n The parameters vector passed to the constraints function\nin the single shooting case must have 6 elements\n');
end

c = [];  % empty inequality constraints.

% define variables for better reading
x1 = state1(1);
y1 = state1(2);
Vx1 = state1(3);
Vy1 = state1(4);
t_i = state1(5);
t_f = state1(6);
mu = parameters(1);        % mass ratio
r_E = parameters(5);       % height of Earth's parking orbit
r_M = parameters(6);       % height of Moon's target orbit

% compute the state at final time:
if nargout < 3
    state2 = PBRFBP_flow(t_i,t_f,state1(1:4),parameters(1:4));  % only state at final time required
else
    [state2,STM] = PBRFBP_STM_variational(t_i,t_f,state1(1:4),parameters(1:4)); % also the STM to be computed
end
x2  = state2(1);
y2  = state2(2);
Vx2 = state2(3);
Vy2 = state2(4);

c_eq = [ (x1 + mu)^2 + y1^2 - r_E^2;                                % distance from Earth at t_i
         (x1 + mu)*(Vx1 - y1) + y1*(Vy1 + x1 + mu);                 % V0 tangent to Earth's parking orbit
         (x2 + mu - 1)^2 + y2^2 - r_M^2;                          % distance from Moon at t_f
         (x2 + mu - 1)*(Vx2 - y2) + y2*(Vy2 + x2 + mu - 1)];  % Vf tangent to moon's parking orbit

if nargout > 2
    Gc = [];    % inequality constraints are not present
    Gc_eq = zeros(4,6);    % initialize jacobian matrix
    Gc_eq(1:2,1:4) = [ 2*(x1 + mu), 2*y1, 0, 0;
                        Vx1, Vy1 , x1 + mu, y1];   % derivative of Psi_1 with respect to initial state
    Dpsi2 = [ 2*(x2 + mu - 1), 2*y2, 0, 0;
              Vx2, Vy2, x2 + mu - 1, y2];          % derivative of Psi_2 with respect to state at t_f
    Gc_eq(3:4,1:4) = Dpsi2*STM;                    % derivative of Psi2 with respect to original state
    Gc_eq(3:4,5) = -Dpsi2*STM*PBRFBP_ode(t_i,state1(1:4),parameters(1:4));
    Gc_eq(3:4,6) = Dpsi2*PBRFBP_ode(t_f,state2(1:4),parameters(1:4));
    Gc_eq = Gc_eq';
end

end


%% POINT 3 FUNCTIONS


function[F,G] = Earth_Moon_cost_ms(y,parameters)
% Function that computes the total cost associated to the two impulse Earth-Moon transfer specified
% by the values contained in the vector y. The initial and final delta velocity are computed with
% respect to prescribed and fixed circular orbits. Multiple shooting technique adopted.
%
% INPUT:
% y          = 4N + 2 augmented state vector built as: (state1,state2,...,stateN, t_i, t_f). [-]
% parameters = vector of parameters built as: 
%         mu = mass ratio governing the P1-P2 motion       [-]
%         mass_S = scaled mass of P3 (Sun)                 [-]    (not used here)
%         rho = scaled P3-(P1-P2) distance                 [-]    (not used here)
%         omega_S = scaled angular velocity of P3 (Sun)    [-]    (not used here)
%         r_E = height of Earth's parking orbit            [-]
%         r_M = height of Moon's parking orbit             [-]
%         Re  = Earth's radius divided by the length unit  [-]    (not used here)
%         Rm  = Moon's radius divided by length unit       [-]    (not used here)
% Some parameters in the vector won't be used but reported to have consistancy with input 
% of the constraints function.
% OUTPUT:
% F    = scalar cost index associated to the total cost of the transfer defined as
%          |v_i - v1| + |v_f - vN|                               [-]
% G    = Gradient of the cost funtion with respect to the NLP variables. 6x1 vector. [-]
%        Computed only when required.

% check on the size of the parameters vector 
if length(parameters) ~= 8
    error('\n The parameters vector passed to the objective function\nin the multiple shooting case must have 8 elements\n');
end

% retrieve state at initial and final points
x1 = y(1);
y1 = y(2);
Vx1 = y(3);
Vy1 = y(4);
xN = y(end - 5);
yN = y(end - 4);
VxN = y(end - 3);
VyN = y(end - 2);
% retrieve parameters for DeltaV computation
mu = parameters(1);
r_E = parameters(5);
r_M = parameters(6);
% compute the initial and final DeltaVs
DeltaV_i = sqrt( (Vx1 - y1).^2 + (Vy1 + x1 + mu).^2 ) - sqrt( (1 - mu)/r_E );
DeltaV_f = sqrt( (VxN - yN).^2 + (VyN + xN + mu - 1).^2 ) - sqrt(mu/r_M);
% overall cost is then:
F = DeltaV_i + DeltaV_f;


% computing the gradient:
if nargout > 1
    G = zeros(length(y),1);     % initialize gradient
    G(1:4) = 1/sqrt( (Vx1 - y1).^2 + (Vy1 + x1 + mu).^2 ).*[ Vy1 + x1 + mu;
                                                              y1 - Vx1;
                                                              Vx1 - y1;
                                                              Vy1 + x1 + mu];                    % components related to state1
    G(end-5:end-2) = 1/sqrt( (VxN - yN).^2 + (VyN + xN + mu - 1).^2 ).*[ VyN + xN + mu - 1;
                                                                         VyN - VxN;
                                                                         VxN - yN;
                                                                         VyN + xN + mu - 1];     % components related to stateN
end
end



function [C,C_eq,Gc,Gc_eq] = Earth_Moon_nonlincon_ms(y,parameters)
% Function expressing the nonlinear constraints for the Earth-Moon transfer optimization using the
% simple shooting technique. The 4 scalar constraints to be enforced concerns the physical distance
% and the velocity direction at both initial and final time. Multiple shooting method considered.
%
% INPUT:
% state1 = 4N + 2 augmented state vector built as: (state1,state2,...,stateN, t_i, t_f). [-]
% parameters = vector of parameters built as: 
%         mu = mass ratio governing the P1-P2 motion       [-]
%         mass_S = scaled mass of P3 (Sun)                 [-]
%         rho = scaled P3-(P1-P2) distance                 [-]
%         omega_S = scaled angular velocity of P3 (Sun)    [-]
%         r_E = height of Earth's parking orbit            [-]
%         r_M = height of Moon's parking orbit             [-]
%         Re  = Earth's radius divided by the length unit  [-]
%         Rm  = Moon's radius divided by length unit       [-]
% OUTPUT:
% c     = nonlinear inequality constraints vector of length 2N+1                [-] 
% c_eq  = 4N vectorial funtion expressing nonlinear equality constraints.       [-]
% Gc    = Jacobian of nonlinear inequality constraints having size (2N+1,4N+2)  [-]
% Gc_eq = Jacobian of nonlinear equality constraints having size (4N,4N+2)      [-]

% check on the size of the parameters vector 
if length(parameters) ~= 8
    error('\n The parameters vector passed to the constraints function\nin the multiple shooting case must have 8 elements\n');
end

% retrieve number of points and other parameters:
N = (length(y)-2)/4;
t1 = y(end - 1);
tN = y(end);
x1 = y(1:4);                % state 1
xN = y(end - 5:end - 2);    % state N
mu = parameters(1);
r_E = parameters(5);        % Earth's parking orbit
r_M = parameters(6);        % Moon's parking orbit
Re  = parameters(7);        % Earth's radius
Rm  = parameters(8);        % Moon's radius
% define time grid:
t_vec = linspace(t1,tN,N);
% initialize vectors of constraints
C = zeros(2*N + 1,1);
C_eq = zeros(4*N,1);

if nargout > 2     % initialize also Jacobian matrixes
    Gc = zeros(2*N + 1,4*N + 2);
    Gc_eq = zeros(4*N,4*N + 2);
end

% start the loop to build the constraints vectors
for i = 1:N - 1
    t_i = t_vec(i);                             % time i
    t_ip1 = t_vec(i+1);                         % time i+1
    state_i = y(4*i - 3 : 4*i);                 % state i
    state_ip1 = y(4*(i+1) - 3 : 4*(i+1));       % state i+1
    
    if nargout < 3 % Jacobian not requested
        flow_of_i = PBRFBP_flow(t_i,t_ip1,state_i,parameters(1:4));                    % compute only flow from state i to state i+1
    else           % Jacobian requested
        [flow_of_i,STM_i] = PBRFBP_STM_variational(t_i,t_ip1,state_i,parameters(1:4)); % compute both flow and STM from state i to state i+1
        Gc_eq(4*i-3 : 4*i , 4*i-3 : 4*(i+1)) = [STM_i, -eye(4)];                       % corresponding lines of equality constraints Jacobian

        Gc_eq(4*i-3 : 4*i , end-1) = -(N - i)/(N - 1)*STM_i*PBRFBP_ode(t_i,state_i,parameters(1:4)) + ...
                                      (N - i - 1)/(N - 1)*PBRFBP_ode(t_ip1,state_ip1,parameters(1:4));   % derivative of z_i with respect to t1
        Gc_eq(4*i-3 : 4*i , end)   = -(i - 1)/(N - 1)*STM_i*PBRFBP_ode(t_i,state_i,parameters(1:4)) + ...
                                       i/(N - 1)*PBRFBP_ode(t_ip1,state_ip1,parameters(1:4));            % derivative of z_i with respect to tN

        Gc(2*i-1 : 2*i , 4*i-3 : 4*i) = [-2*(state_i(1) + mu), -2*state_i(2), 0, 0;
                                         -2*(state_i(1) + mu - 1), -2*state_i(2), 0, 0]; % corresponding lines of inequality constraints Jacobian
    end
    C(2*i-1 : 2*i) = [ Re^2 - (state_i(1) + mu).^2 - state_i(2).^2;
                       Rm^2 - (state_i(1) + mu).^2 - state_i(2).^2];  % inequality constraint to avoid impact on bodies
    C_eq(4*i-3 : 4*i) = flow_of_i - state_ip1;                        % equality constraint z_i
end

% now fill the entries for vectors ad matrixes corresponding to indexes equal to N
% which are not addressed inside the previous loop
C(end-2:end) = [ Re^2 - (xN(1) + mu).^2 - xN(2).^2;            % Earth impact at node N
                 Rm^2 - (xN(1) + mu - 1).^2 - xN(2).^2;        % Moon impact at node N
                 t1 - tN];                                     % positive time of flight
C_eq(end-3:end) = [ (x1(1) + mu).^2 + x1(2).^2 - r_E^2;                                    % initial position on Earth's orbit
                    (x1(1) + mu).*(x1(3) - x1(2)) + x1(2).*(x1(4) + x1(1) + mu);           % initial velocity along Earth's orbit
                    (xN(1) + mu - 1).^2 + xN(2).^2 - r_M.^2;                                  % final position on Moon's orbit
                    (xN(1) + mu - 1).*(xN(3) - xN(2)) + xN(2).*(xN(4) + xN(1) + mu - 1)];  % final velocity along Moon's orbit
if nargout > 3
    Gc(end-2:end-1 , 4*N-3:4*N) = [-2*(xN(1) + mu), -2*xN(2), 0, 0;
                                   -2*(xN(1) + mu - 1), -2*xN(2), 0, 0];       % avoid impacts for node N
    Gc(end , end-1:end) = [ 1, -1];        % related to the time inequality constraint
    Gc_eq(end-3:end-2, 1:4)       = [2*(x1(1) + mu), 2*x1(2), 0, 0;
                                        x1(3), x1(4), x1(1) + mu, x1(2)];     % derivative of Psi1 w.r.t. state1
    Gc_eq(end-1:end, end-5:end-2) = [2*(xN(1) + mu - 1), 2*xN(2), 0, 0;
                                        xN(3), xN(4), xN(1) + mu - 1, xN(2)]; % derivative of PsiN w.r.t stateN
    Gc = Gc';
    Gc_eq = Gc_eq';   % Matlab requires the Jacobian to be transposed
end
end






