%% ASSIGNMENT 1 EX 1

% initialize kernels and other settings
close all;
clear; clc;

% change plots appearance 
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

%% EX1.1
mass_ratio = 3.0359e-6;
% plot the OMx function with respect to x, having set y=z=0
r1_abs = mass_ratio;       % distance of larger body from origin
r2_abs = 1 - mass_ratio;   % distance of secondary body from origin

L3_x_coord = linspace(-1.5,-r1_abs-0.01,1000);
L1_x_coord = linspace(0.5,r2_abs-1e-05,1000);
L2_x_coord = linspace(1.0001,1.5,1000);
L3_points = Collinear_points(L3_x_coord,mass_ratio);
L1_points = Collinear_points(L1_x_coord,mass_ratio);
L2_points = Collinear_points(L2_x_coord,mass_ratio);

figure(); % plot the whole function along with the two bodies
grid on;
hold on;
plot(L3_x_coord,L3_points, L1_x_coord,L1_points, L2_x_coord,L2_points,'LineWidth',3,'Color',[0.9290 0.6940 0.1250]);
scatter(-r1_abs,0,55,[0.8500 0.3250 0.0980],'o','filled');
scatter(r2_abs,0,55,'ob','filled');
xline(0,'k','LineWidth',0.7); yline(0,'k','LineWidth',0.7);
ylim([-3 3]);
xlabel('$x$ [-]'); ylabel('$U$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('$U\left(x\right)$','','','$m_1$','$m_2$','FontSize',15);

figure(); % plot only second and third branches 
grid on;
hold on;
plot(L1_x_coord,L1_points, L2_x_coord,L2_points,'LineWidth',3,'Color',[0.9290 0.6940 0.1250]);
scatter(r2_abs,0,55,'ob','MarkerFaceColor','b');
xline(1,'k','LineWidth',0.7); yline(0,'k','LineWidth',0.7);
xlim([0.95 1.05]);
ylim([-2 2]);
xlabel('$x$ [-]'); ylabel('$U$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('$U\left(x\right)$','','$m_2$','FontSize',15);

% find the x-coordinate for the L2 point
opt1 = optimset('Display','iter','TolX',1e-14);
[x_L2_fzero, f_val, exitflag, ~] = fzero( @(x) Collinear_points(x,mass_ratio), [1.0016 1.02], opt1);  % use fzero
fprintf('\nThe x coordinate of L2 computed with fzero is: %.13e\n', x_L2_fzero);
%% EX 1.2

% initial conditions:
x0 = [1.008296144180133; 0; 0.001214294450297; 0; 0.010020975499502; 0];

% validate the equations of motions through Jacobi's constant
r1_fun = @(x) sqrt( (x(1) + mass_ratio).^2 + x(2).^2 + x(3).^2 );
r2_fun = @(x) sqrt( (x(1) + mass_ratio - 1).^2 + x(2).^2 + x(3).^2 );
U = @(x) 1/2*( x(1).^2 + x(2).^2) + (1 - mass_ratio)./r1_fun(x) + mass_ratio./r2_fun(x) + 1/2*mass_ratio*(1 - mass_ratio); % potential energy
C = @(x) 2*U(x) - ( x(4).^2 + x(5).^2 + x(6).^2 );   % Jacobi's constant
% propagate initial state for arbitrary time:
opt_validation = odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
[tt_val,xx_val] = ode113( @(t,x) RCTBP_ode(t,x,mass_ratio), [0 3/2*pi], x0, opt_validation);
% compute Jacobi constant for each point
C_vec = zeros(length(xx_val),1);
for i = 1:length(xx_val)
    C_vec(i) = C(xx_val(i,:));
end
figure();
semilogy(tt_val,C_vec(:) - C(x0)*ones(length(tt_val),1),'o','LineWidth',2);
hold on;
grid on;
xlabel('Elapsed time [-]'); ylabel('$\Delta$C');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;

% implement the differential corrector: double loop correcting first v_y to get Vx_e = 0 and then
% Vz_e = 0;
iter = 0;  % number of iterations required
Vx_e = 1;
Vz_e = 1;  % set to one in order to enter the while loop
guess_state0 = [x0; reshape( eye(6), [36 1])]; % set initial condition, we are integrating the 42 variational equations
opt3 = odeset('AbsTol',2.5e-14,'RelTol',2.5e-14,'Events',@y_crossing_event);

fprintf('\nStarting the differential corrector loop of point 1.2');
tic;
while abs(Vz_e) > 1e-7 || abs(Vx_e) > 1e-7

    [~,~,t_e,xx_e,~] = ode113(@(t,x) RCTBP_ode_variational(t,x,mass_ratio), [0 2*pi], guess_state0, opt3); % first integration until event occurs
    STM_e = reshape(xx_e(7:end), [6 6]);               % retrieve STM at event
    Vx_e = xx_e(4);                                    % retrieve Vx at event
    Vz_e = xx_e(6);                                    % retrieve Vz at event
    % set the linear system to be solved
    A = [STM_e(4,1), STM_e(4,5);
         STM_e(6,1), STM_e(6,5)];                      % matrix of coefficients
    b = [-Vx_e ; -Vz_e];                               % known term
    Delta_vec = A\b;                                   % solve the linear system to find (Delta_x, Delta_Vy)
    guess_state0(1) = guess_state0(1) + Delta_vec(1);              % update initial x coordinate
    guess_state0(5) = guess_state0(5) + Delta_vec(2);              % update initial Vy
    iter = iter + 1;                                   % update iteration counter
end

fprintf('\nComputations required %f seconds\n',toc);
% initial state generating the periodic orbit is then:
x0_halo = guess_state0(1:6);
% propagate solution and plot:
[~,xx_halo0] = ode113(@(t,x) RCTBP_ode(t,x,mass_ratio), [0 t_e], x0_halo, opt3);   % propagation required only up to half period
figure();
plot3(xx_halo0(:,1), xx_halo0(:,2), xx_halo0(:,3),'r','Linewidth',5);             % plot the integrated leg
grid on;
hold on;
plot3(xx_halo0(:,1), -xx_halo0(:,2), xx_halo0(:,3),'b','Linewidth',5);            % other leg comes by symmetry on xz plane
plot3(x_L2_fzero,0,0,'ok','MarkerSize',10,'MarkerFaceColor','k');      % plot also L2 point          
legend('First branch','Second branch','$L_2$','FontSize',15,'Location','northwest');
xlabel('$x [-]$'), ylabel('$y [-]$'); zlabel('$z [-]$');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15; ax.ZAxis.FontSize = 15;
%% EX 1.3

% generate a set of z0 coordinates:
N = 15;
z0_vec = linspace(x0(3), 0.0046, N);
% initialize the vector of initial conditions that lead to halo orbits
halo_vec0 = zeros(6,N);             % vector containing initial conditions for the halo orbits
halo_vec0(:,1) = x0_halo;           % first initial condition from previous solution
halo_te = zeros(N,1);               % initialize vector for the times of y crossing (half orbital period)

% start the for loop where z0 is increased at each iteration
fprintf('\nStart computing the family of halo orbits using %d points for z0',N);
tic;
for i = 1:N

    Vz_e = 1;
    Vx_e = 1; % set to one in order to enter the while loop at each value of i
    guess_state = [halo_vec0(:,i); reshape( eye(6), [36 1])];     % initial condition for the integration of variational equation
    count = 1;
    while abs(Vz_e) > 1e-7 || abs(Vx_e) > 1e-7                % start the differential correction procedure for the current z0 value

        [t_halo,x_halo,t_e,xx_e,~] = ode113(@(t,x) RCTBP_ode_variational(t,x,mass_ratio), [0 2*pi], guess_state, opt3); % first integration until event occurs
        STM_e = reshape(xx_e(7:end), [6 6]);               % retrieve STM at event
        Vx_e = xx_e(4);                                    % retrieve Vx at event
        Vz_e = xx_e(6);                                    % retrieve Vz at event
        % set the linear system to be solved
        A = [STM_e(4,1), STM_e(4,5);
             STM_e(6,1), STM_e(6,5)];                      % matrix of coefficients
        b = [-Vx_e ; -Vz_e];                               % known term
        Delta_vec = A\b;                                   % solve the linear system to find (Delta_x, Delta_Vy)
        guess_state(1) = guess_state(1) + Delta_vec(1);              % update initial x coordinate
        guess_state(5) = guess_state(5) + Delta_vec(2);              % update initial Vy
        count = count + 1;
    end

    halo_vec0([1 5],i)   = guess_state([1 5]);             % store results of current iteration               
    halo_te(i) = t_e;                                      % store the time at which the trajectory crosses the xz plane (y_event = 0)
    if i < N
            halo_vec0([1 5],i+1) = guess_state([1 5]);             % update the x0 and Vy0 coordinates for the next iteration as the ones found during this step
            halo_vec0(3,i+1) = z0_vec(i+1);                        % increase the z0 coordinate for the next iteration 
    end
end
fprintf('\nComputations required %f seconds\n',toc);

% try and plot the orbits
f1 = figure(); ax1 = gca; grid on; hold on;
f2 = figure(); ax2 = gca; grid on; hold on;
f3 = figure(); ax3 = gca; grid on; hold on;

opt4 = odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
for i = 1:N
   [~,xx_halo] = ode113(@(t,x) RCTBP_ode(t,x,mass_ratio), [0 halo_te(i)], halo_vec0(:,i), opt4);   % propagation required only up to half period
   position_halo = [xx_halo(:,1:3); flip(xx_halo(:,1)), flip(-xx_halo(:,2)), flip(xx_halo(:,3))];     % other half branch obtained by symmetry
   h1 = plot3(ax1,position_halo(:,1), position_halo(:,2), position_halo(:,3),'Linewidth',5);
   h2 = plot(ax2,position_halo(:,2), position_halo(:,3),'Linewidth',5);
   h3 = plot(ax3,position_halo(:,1), position_halo(:,3),'Linewidth',5);
   if i == 1
       label = 'First orbit';
       h1.DisplayName = label; lab1(1) = h1; 
       h2.DisplayName = label; lab2(1) = h2; 
       h3.DisplayName = label; lab3(1) = h3; 
   elseif i == N
       label = 'Last orbit';
       h1.DisplayName = label; h1.Color = 'r'; lab1(2) = h1; 
       h2.DisplayName = label; h2.Color = 'r'; lab2(2) = h2; 
       h3.DisplayName = label; h3.Color = 'r'; lab3(2) = h3; 
       legend(lab1(:));
       legend(lab2(:));
       legend(lab3(:));
   end
end
legend(ax1,'FontSize',15);
xlabel(ax1,'$x [-]$'); ylabel(ax1,'$y [-]$'); zlabel(ax1,'$z [-]$');
ax1.XAxis.FontSize = 15; ax1.YAxis.FontSize = 15; ax1.ZAxis.FontSize = 15;

legend(ax2,'FontSize',15);
xlabel(ax2,'$y [-]$'); ylabel(ax2,'$z [-]$'); 
ax2.XAxis.FontSize = 15; ax2.YAxis.FontSize = 15;

legend(ax3,'FontSize',15);
xlabel(ax3,'$x [-]$'); ylabel(ax3,'$z [-]$');
ax3.XAxis.FontSize = 15; ax3.YAxis.FontSize = 15;

%% EX1.1 functions
function OMx = Collinear_points(x,mass_ratio)
% function that rules the position of the collinear points for a nondimensional 3D RCTBP along the rotating x axis.
% y and z coordinates are enforced to have 0 value so it becomes a 1D problem.
% INPUT
% x = scalar nondimensional position along x axis [-]
% mass_ratio = parameter ruling the nondimensional RCTBP [-]
%
% OUTPUT
% OMx = (or vx_dot) derivative of the potential function with respect to x [-]
mu = mass_ratio;
OMx = x - (1-mu)./(abs(x + mu).^3).*(x + mu ) - mu./(abs(x + mu - 1).^3).*(x + mu - 1);  % equation for OMx (or vx_dot) having as variable the position x
end

%% EX1.2 functions
function [x_dot] = RCTBP_ode(~,state,mu)
% Function that expresses the ode for the 3D RCTBP
% INPUT:
% ~(t)       = time is not necessary because the system is autonomous [-]
% state      = 6x1 vector containing the current state (r;v)  [-]
% mu         = adimensional paramenter ruling the dynamics    [-]
%
% OUTPUT:
% x_dot  = 6x1 vector for the state derivative (v;v_dot)   [-]

% retrieve parameters
x = state(1);
y = state(2);
z = state(3);
r1 = sqrt( (x + mu).^2 + y.^2 + z.^2 );
r2 = sqrt( (x + mu - 1).^2 + y.^2 + z.^2);

% initialize x_dot(v;acc)
x_dot = zeros(6,1);
% velocity is derivative of position
x_dot(1:3) = state(4:6);
% acceleration terms
x_dot(4) = 2*x_dot(2) + x - (1 - mu)*(mu + x)./(r1.^3) + mu*(1 - mu - x)./(r2.^3);
x_dot(5) = -2*x_dot(1) + y - y*(1 - mu)./(r1.^3) - mu*y./(r2.^3);
x_dot(6) = -(1 - mu)*z./(r1.^3) - mu*z./(r2.^3);
end


function [dXdT] = RCTBP_ode_variational(~,state,mu)
% Function that expresses the ode for the 3D RCTBP and its Jacobiam matrix when requested in order
% to obtain the flow at a given time and the corresponding STM
% INPUT:
% ~(t)       = time is not necessary because the system is autonomous  [-]
% state      = 42x1 vector containing the current state (r;v;STM components) [-]
% mu         = adimensional paramenter ruling the dynamics  [-]
%
% OUTPUT:
% dXdT  = 42x1 vector for the extended state derivative (v;v_dot;STM_dot)  [-]

% retrieve parameters

x = state(1);
y = state(2);
z = state(3);
r1 = sqrt( (x + mu).^2 + y.^2 + z.^2 );
r2 = sqrt( (x + mu - 1).^2 + y.^2 + z.^2);
PHI = reshape(state(7:end), [6 6]);

% initialize x_dot(v;acc)
x_dot = zeros(6,1);
% velocity is derivative of position
x_dot(1:3) = state(4:6);
% acceleration terms
x_dot(4) = 2*x_dot(2) + x - (1 - mu)*(mu + x)./(r1.^3) + mu*(1 - mu - x)./(r2.^3);
x_dot(5) = -2*x_dot(1) + y - y*(1 - mu)./(r1.^3) - mu*y./(r2.^3);
x_dot(6) = -(1 - mu)*z./(r1.^3) - mu*z./(r2.^3);

% we compute also the Jacobian of the dynamics equation. The non trivial terms
% are:
OMxx = 1 - (1 - mu)/(r1.^3) + 3*(1 - mu)*((x + mu).^2)./(r1.^5) - mu/(r2.^3) + 3*mu*((x + mu - 1).^2)./(r2.^5);
OMyy = 1 - (1 - mu)/(r1.^3) + 3*(1 - mu)*(y.^2)./(r1.^5) - mu/(r2.^3) + 3*mu*(y.^2)./(r2.^5);
OMzz = -(1 - mu)/(r1.^3) + 3*(1 - mu).*(z.^2)./(r1^5) - mu/(r2.^3) + 3*mu*(z.^2)./(r2.^5);
OMxy = 3*(1 - mu)*(x + mu).*y./(r1.^5) + 3*mu*(x + mu - 1).*y./(r2.^5);   % equal to OMyx
OMxz = 3*(1 - mu)*(x + mu).*z./(r1.^5) + 3*mu*(x + mu - 1).*z./(r2.^5);   % equal to OMzx
OMyz = 3*(1 - mu).*y.*z./(r1.^5) + 3*mu*y.*z./(r2.^5);                    % equal to OMzx
% the jacobian matrix can now be assembled
A = [zeros(3),eye(3);
     OMxx, OMxy, OMxz, 0, 2, 0;
     OMxy, OMyy, OMyz, -2, 0, 0;
     OMxz, OMyz, OMzz, 0, 0, 0];  % 6x6 matrix

% derivative of the STM:
PHI_dot = A*PHI;
% get back to a vectorial expression of the augmented state
dXdT = [x_dot; PHI_dot(:)];
end

function [x_f] = RCTBP_flow(t_i, t_f, x0, mu)
% Function that computes the flow of the ODE
% INPUT:
% t_i         = initial time of integration      [-]
% t_f         = final time of integration   [-]
% x0          = 6x1 vector containing the initial state state (r0;v0) [-]
% mu          = adimensional paramenter ruling the dynamics   [-]
%
% OUTPUT:
% x_f  = flow of the ode at the final time (x_f, v_f)  [-]
%
% External functions called: RCTBP_ode

% computation of the flow at t_f:
opt = odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
[~,xx] = ode113(@(t,x) RCTBP_ode(t,x,mu), [t_i, t_f], x0, opt);
x_f = xx(end,:)';
end


function [value,isterminal,direction] = y_crossing_event(~,x)
% event function to locate at which time the orbit reaches agian y(t_e) = 0
value = x(2);   % event occurs when y = 0
isterminal = 1; % integration to be stopped when event is met
direction = -1; % event to be located when y is decreasing
end