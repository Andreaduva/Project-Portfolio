%% GROUP 18, EXERCISE 1

clc,clear,close all

set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultaxesticklabelinterpreter','latex')
set(groot,'defaultlegendinterpreter','latex')
% initialize the parallel cluster
pp = gcp;
if pp.Connected ~= 1
   parpool;
end
% INITIALIZE DATA AND VARIABLES

a = 4;  % D'Uva
b = 14; % Nunu
c = 19; % Stabile

omega = 5; % Hz
beta = 0.2; % Hz
v0 = 2; % V

par.v = @(t) v0*cos(omega*t)*exp(-beta*t);
par.k_m = 100/a; % TBD
par.R = 200 + b + c; % Ohm
par.L = 2e-3; % H
par.J1 = 0.5; % Kg*m^2
par.J2 = 0.3; % Kg*m^2
par.b = 0.1; % N/rad/s
par.k = 0.5; % N/rad
t0 = 0; % s
tf = 30; % s
tspan = [t0,tf];
y0 = zeros(5,1);
% SIMULATION USING GIVEN VALUES
[tt,yy] = ode15s(@ ex1_system,tspan,y0,odeset,par);

figure()
plot(tt,yy(:,[2 4]), 'LineWidth',2), grid on
legend('$\theta_1$','$\theta_2$','FontSize',20);
xlabel('$t \left[ s \right]$', 'FontSize',20); ylabel('$ \theta \left[ rad\right]$', 'FontSize', 20);
title('Simulation with initial values', 'FontSize',25);
ax = gca; ax.FontSize = 20;

figure()
plot(tt,yy(:,[3 5]), 'LineWidth',2), grid on
legend('$\dot{\theta}_1$','$\dot{\theta}_2$', 'FontSize', 20);
xlabel('$t \left[ s \right]$', 'FontSize',25); ylabel('$ \dot{\theta} \left[ \frac{rad}{s}\right]$', 'FontSize', 20);
title('Simulation with initial values', 'FontSize',25);
ax = gca; ax.FontSize = 20;

figure()
plot(tt,yy(:,1), 'LineWidth',2), grid on
legend('$i(t)$', 'FontSize',20);
xlabel('$t \left[ s \right]$', 'FontSize',25); ylabel('$ i \left[ A \right]$', 'FontSize', 20);
title('Simulation with initial values', 'FontSize',25);
ax = gca; ax.FontSize = 20;

% STATE MATRIX EIGENVALUES && FREE RESPONSE

sysMatrix = [   -par.R/par.L   0   -par.k_m/par.L   0      0
                   0           0       1            0      0 
              ([par.k_m     -par.k  -par.b         par.k  par.b])/par.J1
                   0           0       0            0      1
              (   [0         par.k   par.b       -par.k  -par.b])/par.J2 ];

eigVal = eig(sysMatrix);

% plot eigenvalues
figure();
plot(eigVal(1:end),'x','LineWidth',3, 'MarkerSize',15); grid on; hold on;
x = xlim ; y = ylim;
area( [0 x(1)], [y(2) y(2)], y(1), "FaceColor",[0.3010 0.7450 0.9330],'FaceAlpha',0.2);
xline(0,'k','LineWidth',2); yline(0,'k','LineWidth',2);
xlabel('Re $(\lambda)$','FontSize',20); ylabel('Im $(\lambda)$', 'FontSize',20);
legend('$\lambda_i$','Region of analytical stability', 'FontSize',20,'Location','nw');
title('State matrix eigenvalues', 'FontSize',25);
ax = gca; ax.FontSize = 20;


% free response of the system
y0_free = [1, 1, 1, 1, 1];
par_free = par;
par_free.v = @(t) 0;

[tt_free, yy_free] = ode15s(@ ex1_system, tspan, y0_free, odeset, par_free);

% plot free response
figure()
plot(tt_free,yy_free(:,[2 4]), 'LineWidth',4), grid on
legend('$\theta_1$','$\theta_2$','FontSize',20);
xlabel('$t \left[ s \right]$', 'FontSize',20); ylabel('$ \theta \left[ rad\right]$', 'FontSize', 20);
title('Free response with non-homogeneous initial conditions', 'FontSize',25);
ax = gca; ax.FontSize = 20;

figure()
plot(tt_free,yy_free(:,[3 5]), 'LineWidth',4), grid on
legend('$\dot{\theta}_1$','$\dot{\theta}_2$', 'FontSize', 55);
xlabel('$t \left[ s \right]$', 'FontSize',20); ylabel('$ \dot{\theta} \left[ \frac{rad}{s}\right]$', 'FontSize', 20);
title('Free response with non-homogeneous initial conditions', 'FontSize',25);
ax = gca; ax.FontSize = 20;

figure()
plot(tt_free,yy_free(:,1), 'LineWidth',4), grid on
legend('$i(t)$', 'FontSize',20);
xlabel('$t \left[ s \right]$', 'FontSize',20); ylabel('$ i \left[ A \right]$', 'FontSize', 20);
title('Free response with non-homogeneous initial conditions', 'FontSize',25);
ax = gca; ax.FontSize = 20;


% LOAD EXPERIMENTAL SAMPLES
[t_input,dth2_input] = readvars('Samples.txt');
tspan = 0:0.1:30;
kR_initial = [par.k_m; par.R];
err_initial = norm(dth2_input - dth2_fun( kR_initial,tspan,y0,par));

% GRID SEARCH APPROACH
tic
sg_points = [200 200]; % number of points along each search grid dimension, to be reduced is computational time is too high
[k_mVec, RVec, err_matrix] = parallel_SearchGrid(kR_initial, tspan, y0, par, sg_points, dth2_input);
% Minimum found using the searchgrid
[min_column, row_index] = min(err_matrix);
[err_searchg, column] = min(min_column);
row = row_index(column);
kRsol_searchg = [k_mVec(column), RVec(row)];

% PLOTTING GRID SEARCH RESULTS
[X,Y] = meshgrid(k_mVec,RVec);
figure()
surf(X,Y,err_matrix,'edgecolor','none')
hold on
plot3(kRsol_searchg(1), kRsol_searchg(2),0.02,'ow','MarkerSize',5, 'LineWidth',5);
xlabel('$K_m$','FontSize',25); ylabel('R','FontSize',25);
xlim([k_mVec(1),k_mVec(end)]); ylim([RVec(1) RVec(end)]);
c= colorbar; axis vis3d; c.Label.String = 'relative error'; c.Label.FontSize = 20;
text_name = "Minimum:"; text_Km = sprintf('$K_{min}$ = %.2f', kRsol_searchg(1));
text_R = sprintf('$R_{min}$ = %.2f', kRsol_searchg(2)); text_err = sprintf('$err_{min}$ = %f', err_searchg);
str = {text_name, text_Km, text_R, text_err};
annotation('textbox', [0.4, 0.5, 0.1, 0.1],'String',str,'interpreter','latex', 'BackgroundColor','w','FontSize',20);
title('Grid search and its minimum', 'FontSize',25);
ax = gca; ax.FontSize = 20;
view(2);
% generic view of the error function shape 
figure()
surf(X,Y,err_matrix,'edgecolor','none')
hold on;
xlabel('$K_m$','FontSize',25); ylabel('R','FontSize',25);
xlim([k_mVec(1),k_mVec(end)]); ylim([RVec(1) RVec(end)]);
c= colorbar; axis vis3d; c.Label.String = 'relative error'; c.Label.FontSize = 20;
title('Error function plot', 'FontSize',25);
ax = gca; ax.FontSize = 20;


% MINIMIZATION USING FMINCON
fmc0 = kR_initial;
opt = optimoptions('fmincon','disp','iter-detailed','PlotFcn','optimplotfval');
[kRsol_fmin,err_fmin,flag,outfmc,lambda,grad] = fmincon(@(kR) norm(dth2_input - dth2_fun(kR,tspan,y0,par)),fmc0,[],[],[],[],[0,0],[],[],opt);

% GENETIC ALGORITHMS
ga_iter = 10;
population_size = 100;           % parameters to be reduced if computational time is excessive
max_generation = 13;
ga_fun = @(kR) norm(dth2_input - dth2_fun(kR, tspan,y0,par));
[kRga_vect, err_ga_vect] = parallel_ga_search(ga_fun, population_size, max_generation, ga_iter);
[err_ga, index] = min(err_ga_vect);
kRsol_ga = kRga_vect(index,:);

% PLOTS

par.k_m = kRsol_fmin(1);
par.R = kRsol_fmin(2);
[tt2,yy2] = ode15s(@ ex1_system,tspan,y0,odeset,par);

figure()
plot(t_input,dth2_input,'o',tt2,yy2(:,5), 'LineWidth',2, 'MarkerSize',4); grid on;
legend('experimental','simulated', 'FontSize',20);
xlabel('$ t \left[ s \right]$', 'FontSize', 20); ylabel('$ \dot{\theta}_1  \left[ \frac{rad}{s} \right]$', 'FontSize', 20);
title('Simulation with chosen parameters against experimental samples', 'FontSize',20);
ax = gca; ax.FontSize = 20;

figure()
plot(t_input,abs(dth2_input-yy2(:,5)), 'LineWidth',2);
xlabel('$ t \left[ s \right]$', 'FontSize', 20); ylabel('$ \textrm{Error magnitude}  \left[-\right]$', 'FontSize', 20);
legend('Error$\left(t_k\right)$', 'FontSize',20);
title('Error at each time sample', 'FontSize',20);
ax = gca; ax.FontSize = 20;


%% EXERCISE 2

clc,clear,close all
% DATA

a = 4;  % D'Uva
b = 14; % Nunu
c = 19; % Stabile

data.water.rho = 1000; % kg/m^3
data.water.c = 4186; % J/kg/K

% Lines

data.line.kcv = 2;
data.line.D = 20e-3; % m
data.line.LT1 = 0.5; % m
data.line.fT1 = 0.032;
data.line.L34 = 1.5; % m
data.line.f34 = 0.032;
data.line.L56 = 0.2; % m
data.line.f56 = 0.04;
data.line.L78 = 2.5; % m
data.line.f78 = 0.035;
data.line.L91 = 2.5; % m
data.line.f91 = 0.028;
data.line.L11 = 1; % m
data.line.f11 = 0.032;

data.line.A = pi*data.line.D^2/4;

% Tank

data.tank.p = 0.1e6; % Pa

% Pump

data.pump.N = 9;
data.pump.Dp = 1.5e-2; % m
data.pump.dp = 0.7e-2; % m
data.pump.p_nom = 5*101325; % Pa

data.pump.Dk = 1e-2; % m
data.pump.lc = 10e-2; % m
data.pump.theta_max = deg2rad(20); % rad
data.pump.n = 4000/60; % rps
data.pump.mk = 2; % kg
data.pump.F0 = 5; % N
data.pump.rk = 1; % Ns/m
data.pump.dk = 1e-3; % m
data.pump.kp = 2.5;

data.pump.cc_max = data.pump.lc*tan(data.pump.theta_max); % m
data.pump.Ap = pi*data.pump.Dp^2/4; % m^2
data.pump.Ak = pi*data.pump.Dk^2/4; % m^2
data.pump.Ak_pipe = pi*data.pump.dk^2/4; % m^2
data.pump.hk = (data.pump.p_nom*data.pump.Ak - data.pump.F0)/data.pump.cc_max;

% Distributor

data.distributor.kd = 15;
data.distributor.d0 = 10e-3; % m
data.distributor.command = @(t) 0.25*t + 0.5;

% Cooler
data.cooler.D = 20e-3; % m
data.cooler.Qc = 100; % W

% Filter

data.filter.kf = 35;
data.filter.kl = 2.5/100;

% Heat exchanger

data.HeatEx.D = 20e-3; % m
data.HeatEx.Le = 0.5; % m
data.HeatEx.Ae = 1000e-4; % m^2
data.HeatEx.k1 = 395; % W/m/K
data.HeatEx.l1 = 1e-2; % m
data.HeatEx.k2 = 310; % W/m/K
data.HeatEx.rho2 = 8620; % kg/m^3
data.HeatEx.c2 = 10*a; % J/kg/K
data.HeatEx.l2 = 2.5e-2; % m
data.HeatEx.k3 = 125; % W/kg/K
data.HeatEx.l3 = 1e-2; % m
data.HeatEx.h = 20; % W/m^2/K
data.HeatEx.k = 100*c; % W

D3 = data.HeatEx.D + 2*data.HeatEx.l3;
D2 = D3 + 2*data.HeatEx.l2;
D1 = D2 + 2*data.HeatEx.l1;
vol = pi*(D2^2 - D3^2)/4*data.HeatEx.Le;
data.HeatEx.C = data.HeatEx.c2*data.HeatEx.rho2*vol;

data.HeatEx.R1 = log(D1/D2)/(2*pi*data.HeatEx.k1*data.HeatEx.Le);
data.HeatEx.R2 = log(D2/D3)/(2*pi*data.HeatEx.k2*data.HeatEx.Le);
data.HeatEx.R3 = log(D3/data.HeatEx.D)/(2*pi*data.HeatEx.k3*data.HeatEx.Le);
data.HeatEx.Rconv = 1/data.HeatEx.Ae/data.HeatEx.h;

% Temperature source

T0 = 400; % K
kT = 2*b;
omega = 5; % 1/s

data.Tsource = @(t) T0 + kT*cos(omega*t);

% SIMULATION

t0 = 0; % s
tf = 25; % s

tspan = [t0,tf];

y0 = [0,0,340,10e-3]';

[tt,yy] = ode45(@ex2_system, tspan, y0, odeset, data);

parout = zeros(length(tt),19);
for ii = 1:length(tt)
[~,parout(ii,:)] = ex2_system(tt(ii),yy(ii,:),data);
end

% PLOTS

figure();
plot(tt,yy(:,1:2), 'LineWidth',2);
legend('\(x \left[m\right] \)','\(v \left[\frac{m}{s}\right]\)','FontSize',25), xlabel(' \(t \left[s\right] \)','FontSize',25), grid on;
ax = gca; ax.FontSize = 25;

figure();
plot(tt,parout(:,end-2),tt,yy(:,3), 'LineWidth',2);
legend('$T_{source}$','$T_2$','FontSize',25), xlabel('\(t \left[s\right]\)','FontSize',25), ylabel('\(\left[K\right]\)','FontSize',25); grid on;
ax = gca; ax.FontSize = 25;

figure();
plot(tt,yy(:,4), 'LineWidth',2);
legend('Tank volume','FontSize',25), xlabel('\(t \left[s\right]\)','FontSize',25), ylabel('$\left[m^3\right]$','FontSize',25); grid on;
ax = gca; ax.FontSize = 25;

figure();
plot(tt,parout(:,end), 'LineWidth',2);
legend('$\Delta T$','FontSize',55); xlabel('\(t \left[s\right]\)','FontSize',25); ylabel('\(\left[K\right]\)','FontSize',25); grid on;
ax = gca; ax.FontSize = 25;

figure();
plot(tt,parout(:,end-1)-273,tt,parout(:,end-1)+parout(:,end)-273, 'LineWidth',2);
legend('$T_{low}$','$T_{high}$','FontSize',25); xlabel('\(t \left[s\right]\)','FontSize',25); ylabel('$\left[^\circ C\right]$','FontSize',25); grid on;
ax = gca; ax.FontSize = 25;

figure();
plot(tt,parout(:,8),tt,parout(:,9), 'LineWidth',2);
legend('\(P_{in}\)','\(P_{out}\)','FontSize',25); xlabel('\(t \left[s\right]\)','FontSize',25); ylabel('\(Pressure \left[Pa \right] \)','FontSize',25); grid on;
title('Cooler','FontSize',25);
ax = gca; ax.FontSize = 25;

figure();
plot(tt,parout(:,12),tt,parout(:,13), 'LineWidth',2);
legend('\(P_{in} \)','\(P_{out} \)','FontSize',25); xlabel('\(t \left[s\right]\)','FontSize',25); ylabel('\(Pressure \left[Pa \right] \)','FontSize',25); grid on;
title('Heater','FontSize',25);
ax = gca; ax.FontSize = 25;

figure();
plot(tt,parout(:,1), 'LineWidth',2);
legend('Flow rate','FontSize',25); xlabel('\(t \left[s\right]\)','FontSize',25); ylabel('$\left[\frac{m^3}{s}\right]$','FontSize',25); grid on;
ax = gca; ax.FontSize = 25;

% OPTIMIZING FLUID TEMPERATURE

Qc_max = (400-293.15)/(data.HeatEx.R1+data.HeatEx.R2+data.HeatEx.R3+data.HeatEx.Rconv)/(1-data.filter.kl);

B = data.HeatEx.C*(data.HeatEx.R1 + data.HeatEx.R2/2);

t_transient = 5*B;

Amat = [1 1;
 exp(-t_transient/B) 1];
bmat = [100;Qc_max];

xmat = Amat\bmat;

d = xmat(1);
e = xmat(2);

data.cooler.Qc = @(t) (d*exp(-t/B)+e).*(t<=t_transient)+Qc_max.*(t>t_transient);

[tt_opt,yy_opt] = ode45(@ex2_system, tspan, y0, odeset, data);

parout_opt = zeros(length(tt_opt),19);
for ii = 1:length(tt_opt)
[~,parout_opt(ii,:)] = ex2_system(tt_opt(ii),yy_opt(ii,:),data);
end

figure();
plot(tt_opt,data.cooler.Qc(tt_opt), 'LineWidth',2);
legend('$\dot{Q}_{out}$','FontSize',25); xlabel('\(t \left[s\right]\)','FontSize',25); ylabel('$\left[W\right]$','FontSize',25); grid on;
ax = gca; ax.FontSize = 25;

figure();
plot(tt_opt,parout_opt(:,end-1)-273.15,tt_opt,parout_opt(:,end-1)+parout_opt(:,end)-273.15, 'LineWidth',2);
legend('$T_{low}$','$T_{high}$','FontSize',25); xlabel('\(t \left[s\right]\)','FontSize',25); ylabel('$\left[^\circ C\right]$','FontSize',25); grid on;
ax = gca; ax.FontSize = 25;
title('Fluid temperature with variable $Q_{out}$, $l_2 = 2.5 \, cm$','FontSize',25)

% SECOND ITERATION: THICKER HEAT EXCHANGER

data.HeatEx.l2 = 2e-2 + 2.5e-2; % extending the second layer
D2 = D3 + 2*data.HeatEx.l2;
D1 = D2 + 2*data.HeatEx.l1;
vol = pi*(D2^2 - D3^2)/4*data.HeatEx.Le;
data.HeatEx.C = data.HeatEx.c2*data.HeatEx.rho2*vol;

data.HeatEx.R1 = log(D1/D2)/(2*pi*data.HeatEx.k1*data.HeatEx.Le);
data.HeatEx.R2 = log(D2/D3)/(2*pi*data.HeatEx.k2*data.HeatEx.Le);

Qc_max = (400-293.15)/(data.HeatEx.R1+data.HeatEx.R2+data.HeatEx.R3+data.HeatEx.Rconv)/(1-data.filter.kl);

B = data.HeatEx.C*(data.HeatEx.R1 + data.HeatEx.R2/2);

t_transient = 5*B;

Amat = [1 1;
 exp(-t_transient/B) 1];
bmat = [100;Qc_max];

xmat = Amat\bmat;

d = xmat(1);
e = xmat(2);

data.cooler.Qc = @(t) (d*exp(-t/B)+e).*(t<=t_transient)+Qc_max.*(t>t_transient);

[tt_opt,yy_opt] = ode45(@ex2_system, tspan, y0, odeset, data);

parout_opt = zeros(length(tt_opt),19);
for ii = 1:length(tt_opt)
[~,parout_opt(ii,:)] = ex2_system(tt_opt(ii),yy_opt(ii,:),data);
end

figure();
plot(tt_opt,parout_opt(:,end-1)-273.15,tt_opt,parout_opt(:,end-1)+parout_opt(:,end)-273.15, 'LineWidth',2);
legend('$T_{low}$','$T_{high}$','FontSize',25); xlabel('\(t \left[s\right]\)','FontSize',25); ylabel('$\left[^\circ C\right]$','FontSize',25); grid on;
ax = gca; ax.FontSize = 25;
title('Fluid temperature with variable $Q_{out}$, $l_2 = 4.5 \, cm$','FontSize',25);


%% FUNCTIONS EX1 


function [dy] = ex1_system(tt,yy,par)
%state function
I = yy(1);
th1 = yy(2);
dth1 = yy(3);
th2 = yy(4);
dth2 = yy(5);

dy = [(par.v(tt) - par.R*I - par.k_m*dth1)/par.L
      dth1
      (par.k*(th2-th1) + par.b*(dth2 - dth1) + par.k_m*I)/par.J1
      dth2
      (-par.k*(th2-th1) - par.b*(dth2 - dth1))/par.J2];
end



function [dth2] = dth2_fun(kR,tspan,y0,par)
% auxiliary function to be used in fmincon
par.k_m = kR(1);
par.R = kR(2);

[~,yy] = ode15s(@ex1_system,tspan,y0,odeset,par);

dth2 = yy(:,5);

end



function [k_mVec, RVec, err] = parallel_SearchGrid(kRsol_fmin,tspan,y0, par, sg_points, dth2_input)
% function for the searchgrid
D = parallel.pool.DataQueue;
h = waitbar(0, 'Search grid: computing...', 'Units','centimeters');
N = sg_points(1);
p = 1;
afterEach(D, @nUpdateWaitbar );
% parrallel for loop, requires specific syntax
k_mVec = linspace(0,kRsol_fmin(1) +20,sg_points(1));
RVec = linspace(kRsol_fmin(2)-100,kRsol_fmin(2)+100,sg_points(2));
err = zeros(length(RVec),length(k_mVec));
v = par.v; L = par.L; J1 = par.J1; J2 = par.J2; b = par.b; k = par.k;

tic;
parfor ii = 1:length(k_mVec)
    err_row = zeros(1,length(k_mVec));
    par_aux = struct('v',v,'k_m',k_mVec(ii),'L',L,'J1',J1,'J2',J2,'b',b,'k',k);

    for jj = 1:length(RVec)

        par_aux.R = RVec(jj);
        [~,yy] = ode15s(@ ex1_system,tspan,y0,odeset,par_aux);

        err_row(jj) = norm(dth2_input-yy(:,5));
    end
    err(:,ii) = err_row;
    send(D,ii);
end


    function nUpdateWaitbar(~)
        waitbar(p/N, h);
        p = p + 1;
    end
waitbar(1,h,'Computations concluded');
pause(1);
close(h);
fprintf('Search grid took %.2f seconds.\n',toc);
end



function [sol_vect, fval_vect] = parallel_ga_search(ga_fun, population, generations, n_iter)
% function for ga iterations
D = parallel.pool.DataQueue;
h = waitbar(0, 'Genetic algorithm: computing...', 'Units','centimeters');
N = n_iter;
p = 1;
afterEach(D, @nUpdateWaitbar );
tic;
fval_vect = zeros(n_iter,1);
sol_vect = zeros(n_iter,2);
parfor count = 1:n_iter
    ga_opt = optimoptions('ga','UseParallel', false, 'MaxGenerations', generations,'PopulationSize',population, 'Display','off');
    [sol, fval] = ga( ga_fun, 2, [], [], [], [], [0 0], [30 400], [],ga_opt);
    sol_vect(count,:) = sol; fval_vect(count) = fval;
    send(D,count);
end

    function nUpdateWaitbar(~)
        waitbar(p/N, h);
        p = p + 1;
    end
waitbar(1,h,'Computations concluded');
pause(1);
close(h);
fprintf('Ga took %.2f seconds.\n', toc);
end

%% FUNCTIONS EX 2

function [yy,parout] = ex2_system(tt,xx,data)

xk = xx(1);
vk = xx(2);
T2 = xx(3);
% Vtank = xx(4);

if xk < 0
    xk = 0;
end
if xk <= 0 && vk < 0
    xk = 0;
    vk = 0;
end
if xk > data.pump.cc_max
    xk = data.pump.cc_max;
end
if xk >= data.pump.cc_max && vk > 0
    xk = data.pump.cc_max;
    vk = 0;
end

s = data.pump.dp/data.pump.lc*(data.pump.cc_max - xk);

vp = vk*data.pump.Ak/data.pump.Ak_pipe;

z = data.distributor.command(tt);

if z > 1
    z = 1;
elseif z < 0
    z = 0;
end

alpha = 2*acos(1-abs(2*z));

A_dist = data.distributor.d0^2/8*(alpha-sin(alpha));

% Flow rates

Q1 = data.pump.n*data.pump.N*data.pump.Ap*s;
Q9 = Q1*(1-data.filter.kl);

% Temperatures

Tsource = data.Tsource(tt);

if isa(data.cooler.Qc,'double')
    Q_out = data.cooler.Qc; % power exiting the fluid
else
    Q_out = data.cooler.Qc(tt);
end

dTcooler = Q_out/(data.water.c*data.water.rho*Q1);

Q_in = Q_out*(1-data.filter.kl); % power entering the fluid

T2_dot = ((Tsource-T2)/(data.HeatEx.R1 + data.HeatEx.R2/2) - Q_in)/data.HeatEx.C;

Tf = T2 - Q_in*(data.HeatEx.R2/2+data.HeatEx.R3+data.HeatEx.Rconv);

% Pressures

P1 = data.tank.p - 1/2*data.line.fT1*data.water.rho*data.line.LT1/...
    data.line.D*abs(Q1)*Q1/data.line.A^2;

P11 = data.tank.p + 1/2*data.line.f11*data.water.rho*data.line.L11/...
    data.line.D*abs(Q9)*Q9/data.line.A^2;

P10 = P11/exp(Q_in/data.HeatEx.k);

P9 = P10 + 1/2*data.line.f91*data.water.rho*data.line.L91/...
    data.line.D*abs(Q9)*Q9/data.line.A^2;

P8 = P9 + 1/2*data.filter.kf*data.water.rho*abs(Q9)*Q9/data.line.A^2;

P7 = P8 + 1/2*data.line.f78*data.water.rho*data.line.L78/...
    data.line.D*abs(Q1)*Q1/data.line.A^2;

P6 = P7/exp(-Q_out/data.HeatEx.k);

P5 = P6 + 1/2*data.line.f56*data.water.rho*data.line.L56/...
    data.line.D*abs(Q1)*Q1/data.line.A^2;

P4 = P5 + 1/2*data.distributor.kd*data.water.rho*abs(Q1)*Q1/A_dist^2;

P3 = P4 + 1/2*data.line.f34*data.water.rho*data.line.L34/...
    data.line.D*abs(Q1)*Q1/data.line.A^2;

P2 = P3 + 1/2*data.line.kcv*data.water.rho*abs(Q1)*Q1/data.line.A^2;

Pk = P2 - 1/2*data.pump.kp*data.water.rho*abs(vp)*vp;

% Pressure regulator

xk_dot = vk;
vk_dot = 1/data.pump.mk*(Pk*data.pump.Ak - data.pump.F0 - data.pump.hk*xk - data.pump.rk*vk);
Vtank_dot = Q9 - Q1;

% state equation
yy = [xk_dot
      vk_dot
      T2_dot
      Vtank_dot
      ];

parout = [Q1 Q9 ...
          P1 P2 P3 P4 P5 P6 P7 P8 P9 P10 P11 ...
          vk_dot,A_dist,z,Tsource,Tf,dTcooler];

end