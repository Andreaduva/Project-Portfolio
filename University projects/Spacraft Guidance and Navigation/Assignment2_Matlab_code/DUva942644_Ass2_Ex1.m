%% ASSIGNMENT 2 EX 1

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

% load spice kernels from meta-kernel - NB to be run while in the Assignment01 folder as current
% working directory
% cspice_furnsh('..\assignment02_win.tm');  % Windows version
cspice_furnsh('../assignment02_mac.tm');  % MacOs version
kernel_number = cspice_ktotal('ALL');
fprintf('\nNumber of kernels loaded is %d\n',kernel_number);

% set constants and parameters:
x0_ref = [6054.30795817484; -3072.03883303992; -133.115352431876];
v0_ref = [4.64750094824087; 9.18608475681236; -0.62056520749034];
state_ref_per = [x0_ref;v0_ref];
P0_ref = [5.6e-3  3.5e-3  -7.1e-4   0   0   0;
          3.5e-3  9.7e-3   7.6e-4   0   0   0;
          -7.1e-4 7.6e-4   8.1e-4   0   0   0;
          0       0        0        2.8e-7 0 0;
          0       0        0        0  2.7e-7 0;
          0       0        0        0   0   9.6e-8];
t_ref = cspice_str2et('2022-11-11 19:08:49.824 UTC');
mu = cspice_bodvrd('EARTH','GM',1);
r_Earth = cspice_bodvrd('EARTH','RADII',3);
r_Earth = r_Earth(1);
keplerian_state_per = cspice_oscelt(state_ref_per,t_ref,mu);
a_0 = keplerian_state_per(1)/(1 - keplerian_state_per(2));  % semi-major axis
T   = 2*pi*sqrt( (a_0^3)/mu );                % orbital period [s]
i = rad2deg(keplerian_state_per(3));
fprintf('The reference orbit has an inclination of %.3f degrees.\n',i);

% nominal state at apocenter:
keplerian_state_apo = keplerian_state_per;
state_ref_apo = cspice_conics(keplerian_state_apo,t_ref + T/2);

% plot the nominal orbit
opt_ode = odeset('AbsTol',2.5e-14,'RelTol',2.5e-14);
[tt,yy] = ode113(@(t,y) eom_2body(t,y,mu), [0 T], state_ref_per,opt_ode);

% XY plane view
figure();
hold on;
grid on;
[X,Y,Z] = sphere;
X = r_Earth*X; Y = r_Earth*Y; Z = r_Earth*Z;
surf(X,Y,Z,'FaceColor',[0.5 0.5 0.5]);
plot3(yy(:,1),yy(:,2),yy(:,3),'r','LineWidth',3);
axis equal;
xlabel('$x$ [km]'); ylabel('$y$ [km]'); zlabel('$z$ [km]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15; ax.ZAxis.FontSize = 15;
legend('','Keplerian orbit','FontSize',15);
title('Arian 5 upper stage orbit','FontSize',15);
view(0,90);

% initial Monte Carlo points and uncertainty distribution:
N = 3500;
init_samples = mvnrnd(state_ref_per,P0_ref,N)';
figure();
grid on;
hold on;
plot(init_samples(1,:),init_samples(2,:),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
sigma_ellipse = Cov_ellipse2D(3, x0_ref(1:2), P0_ref(1:2,1:2));
plot(sigma_ellipse(1,:),sigma_ellipse(2,:),'--r','LineWidth',3);
plot(yy(:,1),yy(:,2),'b','LineWidth',2.5);
plot(x0_ref(1),x0_ref(2),'ok','MarkerSize',10,'MarkerFaceColor','k');
legend('Samples', '$3\sigma$ ellipse','Reference trajectory', 'Mean state','FontSize',15,'Location','southeast');
xlabel('$x$ [km]'); ylabel('$y$ [km]');
title('Initial condition at time $t_0$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
ylim([ -3072.39416970682 -3071.69013187496]); xlim([6053.95573672571 6054.77711419622]);

%% uncertainty propagation
t_span = T/2:T/2:4*T;
% initialize vectors used to store results:
LinCov_m_vec   = zeros(6,length(t_span));
LinCov_P       = zeros(6,6,length(t_span));
LinCov_Prr_vec = zeros(1,length(t_span));
LinCov_Pvv_vec = zeros(1,length(t_span));
LinCov_timout = zeros(1,length(t_span));

UT_m_vec   = zeros(6,length(t_span));
UT_P       = zeros(6,6,length(t_span));
UT_Prr_vec = zeros(1,length(t_span));
UT_Pvv_vec = zeros(1,length(t_span));
UT_timout = zeros(1,length(t_span));

MC_m_vec       = zeros(6,length(t_span));
MC_P           = zeros(6,6,length(t_span));
MC_Prr_vec     = zeros(1,length(t_span));
MC_Pvv_vec     = zeros(1,length(t_span));
MC_samples_vec = zeros(6,N,length(t_span));   % store MC points in 3D matrix
MC_timout = zeros(1,length(t_span));

% loop over times to propagate statistics:
fprintf('Beginning the uncertainty propagation...\n');
startLoop = tic;
for i = 1:length(t_span)
    % linear propagation
    startIter = tic;
    [LinCov_m_vec(:,i),P_aux] = LinCov_2body(state_ref_per,P0_ref,t_span(i),mu);
    LinCov_P(:,:,i) = P_aux;
    LinCov_Prr_vec(i) = sqrt(trace(P_aux(1:3,1:3)));
    LinCov_Pvv_vec(i) = sqrt(trace(P_aux(4:6,4:6)));
    endIter = toc(startIter);
    LinCov_timout(i) = endIter;
    % UT propagation
    startIter = tic;
    [UT_m_vec(:,i),P_aux] = UT_2body(state_ref_per,P0_ref,[1e-3; 2; 0; t_span(i); mu]);
    UT_P(:,:,i) = P_aux;
    UT_Prr_vec(i) = sqrt(trace(P_aux(1:3,1:3)));
    UT_Pvv_vec(i) = sqrt(trace(P_aux(4:6,4:6)));
    endIter = toc(startIter);
    UT_timout(i) = endIter;
    % Monte Carlo simulation
    startIter = tic;
    [MC_samples_vec(:,:,i),MC_m_vec(:,i),P_aux] = MC_2body(N,state_ref_per,P0_ref,t_span(i),mu);
    MC_P(:,:,i) = P_aux;
    MC_Prr_vec(i) = sqrt(trace(P_aux(1:3,1:3)));
    MC_Pvv_vec(i) = sqrt(trace(P_aux(4:6,4:6)));
    endIter = toc(startIter);
    MC_timout(i) = endIter;
end
fprintf('The uncertainty propagation process took %.3f seconds.\n',toc(startLoop));

% plot CPU time behaviour for each method
figure();
h = semilogy(1:8,LinCov_timout,'o','MarkerSize',10);
h.MarkerFaceColor = h.Color;
hold on;
h = semilogy(1:8,UT_timout,'o','MarkerSize',10);
h.MarkerFaceColor = h.Color;
grid on;
h = semilogy(1:8,MC_timout,'o','MarkerSize',10);
h.MarkerFaceColor = h.Color;
legend('LinCov', 'UT','MC','FontSize',15);
xlabel('Elapsed semiperiods $\left[\frac{T}{2}\right]$'); ylabel('CPU time $\left[s\right]$');
title('Time performance comparison','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
%% EX1.2
% draw the covariance ellipses for each apogee and perigee:
for i = 1:length(t_span)
    figure();
    hold on;
    grid on;
    %ref_state = twobody_flow(state_ref0,0,t_span(i),mu);
    MC_ellipse = Cov_ellipse2D(3, MC_m_vec(1:2,i), MC_P(1:2,1:2,i));
    LinCov_ellipse = Cov_ellipse2D(3,LinCov_m_vec(1:2,i), LinCov_P(1:2,1:2,i));
    UT_ellipse = Cov_ellipse2D(3,UT_m_vec(1:2,i), UT_P(1:2,1:2,i));
    % plot routine
    plot(MC_ellipse(1,:),MC_ellipse(2,:),'--k','LineWidth',3);
    plot(LinCov_ellipse(1,:),LinCov_ellipse(2,:),'g','LineWidth',3);
    plot(UT_ellipse(1,:),UT_ellipse(2,:),'--r','LineWidth',3);
    plot(yy(:,1),yy(:,2),'b','LineWidth',2.5);
    % plot also all samples
    plot(MC_samples_vec(1,:,i),MC_samples_vec(2,:,i),'o','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
    legend('$3\sigma$ MC', '$3\sigma$ LinCov', '$3\sigma$ UT', 'Reference Trajectory', 'MC samples','FontSize',15);
    xlabel('$x$ [km]'); ylabel('$y$ [km]');
    ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
    switch i
        case 1
            pass_n = 1;
            x_lim = [-5.855e4 -5.82e4]; y_lim = [2.92e4 2.99e4];
            legend('Location','southwest');
        case 2
            pass_n = 1;
            x_lim = [3500 8000]; y_lim = [-7000 1000];
            legend('Location','northwest');
        case 3
            pass_n = 2;
            x_lim = [-5.841e4 -5.833e4]; y_lim = [2.9e4 3.04e4];
            legend('Location','southeast');
        case 4
            pass_n = 2;
            x_lim = [2000 10000]; y_lim = [-10000 4000];
            legend('Location','northwest');
        case 5
            pass_n = 3;
            x_lim = [-5.86e4 -5.81e4]; y_lim = [2.85e4 3.1e4];
            legend('Location','northwest');
        case 6
            pass_n = 3;
            x_lim = [-2000 12000]; y_lim = [-1.5e4 1e4];
            legend('Location','northwest');
        case 7
            pass_n = 4;
            x_lim = [-5.88e4 -5.79e4]; y_lim = [2.8e4 3.1e4];
            legend('Location','southeast');
        case 8
            pass_n = 4;
            x_lim = [-4000 14000]; y_lim = [-2e4 1.5e4];
            legend('Location','southeast');
    end
    if mod(i,2) == 1
        str = 'Apogee';
        title([str,' $',num2str(pass_n),'$. Elapsed time is $',num2str(t_span(i)/3600),'$ hours.'],'FontSize',15);
    else
        str = 'Perigee';
        title([str,' $',num2str(pass_n),'$. Elapsed time is $',num2str(t_span(i)/3600),'$ hours.'],'FontSize',15);
    end
    xlim(x_lim); ylim(y_lim); 
end



% display on screen the results
LinCov_m_delta = zeros(6,8);
LinCov_m_delta(:,[1,3,5,7]) = LinCov_m_vec(:,[1,3,5,7]) - state_ref_apo;
LinCov_m_delta(:,[2,4,6,8]) = LinCov_m_vec(:,[2,4,6,8]) - state_ref_per;
RowName = {'Apogee  passage 1';'Perigee passage 1';'Apogee  passage 2';'Perigee passage 2';...
         'Apogee  passage 3';'Perigee passage 3';'Apogee  passage 4';'Perigee passage 4'};
ColumnName_LinCov = {'Delta x LinCov','Delta y LinCov','Delta z LinCov','Delta Vx LinCov',...
                     'Delta Vy LinCov','Delta Vz LinCov','Prr LinCov','Pvv LinCov'};
T_LinCov = table(LinCov_m_delta(1,:)',LinCov_m_delta(2,:)',LinCov_m_delta(3,:)',...
                 LinCov_m_delta(4,:)', LinCov_m_delta(5,:)', LinCov_m_delta(6,:)',...
                 LinCov_Prr_vec',LinCov_Pvv_vec',...
                 'RowNames',RowName,'VariableNames',ColumnName_LinCov);
disp(T_LinCov);
% compute norm of delta vectors:
LinCov_delta_r = zeros(1,8);
LinCov_delta_v = zeros(1,8);
for i = 1:length(t_span)
    LinCov_delta_r(i) = norm(LinCov_m_delta(1:3,i));
    LinCov_delta_v(i) = norm(LinCov_m_delta(4:6,i));
end


UT_m_delta = zeros(6,8);
UT_m_delta(:,[1,3,5,7]) = UT_m_vec(:,[1,3,5,7]) - state_ref_apo;
UT_m_delta(:,[2,4,6,8]) = UT_m_vec(:,[2,4,6,8]) - state_ref_per;
ColumnName_UT = {'Delta x UT','Delta y UT','Delta z UT','Delta Vx UT','Delta Vy UT',...
                 'Delta Vz UT','Prr UT','Pvv UT'};
T_UT= table(UT_m_delta(1,:)',UT_m_delta(2,:)',UT_m_delta(3,:)',UT_m_delta(4,:)', ...
            UT_m_delta(5,:)',UT_m_delta(6,:)',UT_Prr_vec',UT_Pvv_vec',...
           'RowNames',RowName,'VariableNames',ColumnName_UT);
disp(T_UT);
% compute norm of delta vectors:
UT_delta_r = zeros(1,8);
UT_delta_v = zeros(1,8);
for i = 1:length(t_span)
    UT_delta_r(i) = norm(UT_m_delta(1:3,i));
    UT_delta_v(i) = norm(UT_m_delta(4:6,i));
end

MC_m_delta = zeros(6,8);
MC_m_delta(:,[1,3,5,7]) = MC_m_vec(:,[1,3,5,7]) - state_ref_apo;
MC_m_delta(:,[2,4,6,8]) = MC_m_vec(:,[2,4,6,8]) - state_ref_per;
ColumnName_MC = {'Delta x MC','Delta y MC','Delta z MC','Delta Vx MC','Delta Vy MC',...
                 'Delta Vz MC','Prr MC','Pvv MC'};
T_MC = table(MC_m_delta(1,:)',MC_m_delta(2,:)',MC_m_delta(3,:)',MC_m_delta(4,:)', ...
             MC_m_delta(5,:)',MC_m_delta(6,:)',MC_Prr_vec',MC_Pvv_vec',...
            'RowNames',RowName,'VariableNames',ColumnName_MC);
disp(T_MC);
% compute norm of delta vectors:
MC_delta_r = zeros(1,8);
MC_delta_v = zeros(1,8);
for i = 1:length(t_span)
    MC_delta_r(i) = norm(MC_m_delta(1:3,i));
    MC_delta_v(i) = norm(MC_m_delta(4:6,i));
end

%% EX1.1 FUNCTIONS

function [x_dot] = eom_2body(~, x, mu)
% eom for restricted 2-body problem without perturbations
% INPUT:
%
% x   = state vector (r, v)                               [km],[km/s]
% t   = final epoch - negligible since autonomous system  [s]
% mu  = gravitational parameter                           [km^2/s^2]
%
% OUTPUT:
%
% x_dot = equation of motion for x_dot
%
r_norm = norm( [x(1) x(2) x(3)]);   % magnitude of position vector
x_dot = [ x(4);
          x(5);
          x(6);
          -mu/(r_norm^3)*x(1);
          -mu/(r_norm^3)*x(2);
          -mu/(r_norm^3)*x(3);];
end


function [dx_dphi] = eom_2body_STM(~, x, mu)
% eom ruling the State Transition matrix in restricted 2-body problem without perturbations
% INPUT:
%
% x   = state vector of the system. 42 components with first 6 being (r,v)
% t   = final epoch - negligible since autonomous system  [s]
% mu  = gravitational parameter                           [km^2/s^2]
%
% OUTPUT:
%
% x_dot = equation of motion for x_dot
%

% start with the six scalar equations related to state dynamics
r_norm = norm( [x(1) x(2) x(3)]);   % magnitude of position vector
x_dot = [ x(4);
          x(5);
          x(6);
          -mu/(r_norm^3)*x(1);
          -mu/(r_norm^3)*x(2);
          -mu/(r_norm^3)*x(3);];

% definition of matrix A. Computed analytically in this simple case
A11 = zeros(3); A12 = eye(3); A22 = zeros(3); % base blocks
A21 = 3*mu/(r_norm^5)*[ x(1)^2, x(1)*x(2), x(1)*x(3); x(1)*x(2), x(2)^2, x(2)*x(3); x(1)*x(3), x(2)*x(3), x(3)^2] - mu/(r_norm^3)*eye(3);
A = [A11, A12; A21, A22];

phi = reshape( x(7:end), [6 6]);   % reshape the state vector x(7:42) into a 6x6 matrix
phi_dot = A*phi;                   % multiply by A                      
dx_dphi = [x_dot; phi_dot(:)];     % retrieve column vector
end

function [x_f, Phi_f] = twobody_flow(x_0,t_0,t_f,mu)
% function that defines the flow of the restricted 2 body model and gives
% back the final state and the State Transition Matrix
%
% INPUT 
% x_0   = initial state (r_0, v_0)   [km],[km/s]
% t_0   = initial epoch              [s]
% t_f   = final epoch                [s]
% mu    = gravitational parameter    [km^3/s^2]
%
% OUTPUT
% x_f     = state at final epoch (r_f, v_f)   [km], [km/s]
% Phi_f   = State transition matrix

options = odeset('reltol', 2.5e-14, 'abstol', 2.5e-14);
switch nargout
    case 1
        [~,x_sol] = ode113(@(t,x) eom_2body(t,x,mu),[t_0, t_f],x_0,options);
        x_f = x_sol(end,:)';
    case 2
        [~,x_sol] = ode113(@(t,x) eom_2body_STM(t,x,mu),[t_0, t_f], [x_0; reshape(eye(6), [36, 1])],options);
        x_f = x_sol(end,1:6)';
        Phi_f = reshape( x_sol(end,7:end)', [6 6]);
end

end

function [mean,P] = LinCov_2body(mean0,P0,tof,mu)
% function that propagates uncertainties along a Keplerian orbit using the linearized dynamics
%
% INPUT:
% mean0  = initial mean of the state, [nx1] vector.
% P0     = initial covariance matrix of the state [n x n] matrix.
% tof    = time for which the transformation is to be propagated
% mu     = gravitational parameter

% OUTPUT:
% sample_mean = estimated mean at the end of the transformation, [nx1] vector.
% sample_P    = estimated covariance, [n x n] matrix

% compute the dynamics flow and related STM:
[mean,STM] = twobody_flow(mean0, 0, tof, mu);
P = STM*P0*STM';
end



function [sample_mean,sample_P] = UT_2body(mean0,P0,parameters)
% function that propagates uncertainties along a Keplerian orbit using the most general 
% form of the Unscented Transform.
%
% INPUT:
% mean0 = initial mean of the state, [nx1] vector.
% P0    = initial covariance matrix of the state [n x n] matrix.
% parameters = set of useful parameter required for the propagation:
%              aplha  = parameter related to spreading of sigma points
%              beta   = parameter delivering a-priori knowledge of the distribution
%              k      = parameter usually set to zero
%              tof    = time for which the transformation is to be propagated
%              mu     = gravitational parameter
%
% OUTPUT:
% sample_mean = estimated mean at the end of the transformation, [nx1] vector.
% sample_P    = estimated covariance, [n x n] matrix

alpha = parameters(1);
beta = parameters(2);
k = parameters(3);
tof = parameters(4);
mu = parameters(5);
n = length(mean0);
lambda = alpha^2*(n + k) - n;

% Step 1: create sigma points
sqrtP = sqrtm((n + lambda)*P0);
chi_points = [mean0, mean0 + sqrtP, mean0 - sqrtP];

% compute weights
W0_m = lambda/(n + lambda);
W0_c = lambda/(n + lambda) + (1 - alpha^2 + beta);
Wi_m = 1/(2*(n + lambda));
Wi_c = Wi_m;
W_vec_m = [W0_m, ones(1,2*n)*Wi_m];
W_vec_c = [W0_c, ones(1,2*n)*Wi_c];
% Step 2: propagate the sigma points
Y = zeros(n, 2*n + 1);
for i = 1:2*n + 1
    Y(:,i) = twobody_flow( chi_points(:,i), 0, tof, mu);
end

% Step 3: compute sample mean and covariance
sample_mean = sum(W_vec_m.*Y,2);
sample_P = W_vec_c.*(Y - sample_mean)*(Y - sample_mean)';
end


%% EX1.2 FUNCTIONS

function [MC_samples,MC_mean,MC_P] = MC_2body(N,mean0,P0,tof,mu)
% function that uses a Monte Carlo method to propagate statistics along a keplerian dynamics.
% It is assumed that the initial probability distribuition is Gaussian
%
% INPUT:
% N = number of samples to be taken from initial distribution.
% mean0 = initial mean of the state, [nx1] vector.
% P0    = initial covariance matrix of the state [n x n] matrix.
% tof = time for which dynamics has to be propagated
% mu = gravitational parameter
%
% OUTPUT:
% MC_samples = [6xN] matrix containing all the propagated samples.
% MC_mean    = final sample mean.
% MC_P       = final sample covariance.

% generate N samples
samples = mvnrnd(mean0,P0,N)';
MC_samples = zeros(length(mean0),N);
% propagate them:
parfor i = 1:N
    MC_samples(:,i) = twobody_flow( samples(:,i), 0, tof, mu);
end

% compute sample mean and covariance from them:
MC_mean = sum(MC_samples,2)/N;
MC_P = (MC_samples - MC_mean)*(MC_samples - MC_mean)'/(N-1);
end


function [ellipse] = Cov_ellipse2D(N_std,mean,P)
% function that computes the 2D covariance ellipse with a given confidence interval
%
% INPUT:
% P = [2x2] covariance matrix
% mean  = mean value around which covariance ellipse is defined
% N_std = number of standard deviations which defines the confidence interval
% OUTPUT:
% ellipse = resulting covariance ellipse

confidence = 2*normcdf(N_std) - 1;       % 
scale = chi2inv(confidence,2);           % chi value for given confidence and degrees of freedom
P_scaled = P*scale;                      % scaled covariance matrix
[V, D] = eig(P_scaled);                  % eigenvalues of scaled covariance matrix
[D, order] = sort(diag(D),'descend');    % sort matrix eigenvalues
D = diag(D);
V = V(:,order);
V_scaled = V*sqrt(D);

% now we can get the ellipse points:
t = linspace(0,2*pi,400);
circle = [cos(t); sin(t)];
ellipse = bsxfun(@plus, V_scaled*circle, mean);
end



