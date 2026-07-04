%% ASSIGNMENT 2 EX 3

% initialize kernels and other settings
close all;
clear; clc;
cspice_kclear();

% change plots appearance 
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% load spice kernels from meta-kernel - NB to be run while in the Assignment01 folder as current
% working directory
% cspice_furnsh('..\assignment02_win.tm');  % Windows version
cspice_furnsh('../assignment02_mac.tm');  % MacOs version
kernel_number = cspice_ktotal('ALL');
fprintf('\nNumber of kernels loaded is %d\n',kernel_number);

%% Set Constants
l = 10;      % s/c length [m].
h = 5;       % s/c heigth [m].
d = 3;       % s/c depth  [m].
dim = [l;h;d];
m = 213000;  % s/c mass [kg].
t0_ref = cspice_str2et('2023-04-01 14:55:12.023 UTC'); % reference epoch
q0 = [0.674156352338764;
      0.223585877389611;
      0.465489474399161;
      0.528055032413102];    % initial quaternion
om0 = [-0.001262427155865;
        0.001204540074343;
       -0.000039180139156];  % initial angular velocity
J = m/12*diag([d^2 + h^2; l^2 + h^2; d^2 + l^2]);
a = 42164;                   % Geostationary orbit semi-major axis [km]
mu = cspice_bodvrd('EARTH','GM',1);
n = sqrt(mu/(a^3));          % target spacecraft's mean motion [1/sec]
t_span = t0_ref:1:t0_ref+cspice_spd; % one day span from t0 with samples after each second.
tof_span = t_span - t0_ref;
T = 2*pi/n;
n_periods = 3;
t_span_extended = t0_ref:1:t0_ref+n_periods*T;
tof_span_extended = t_span_extended-t0_ref;
% integrate the kinematics
opt = odeset('AbsTol',2.5e-14,'RelTol',2.3e-14);
[~,yy] = ode113( @(t,y) kinematics_eom(t,y,J), tof_span, [om0;q0],opt);
om_vec = yy(:,1:3)';
q_vec  = yy(:,4:7)';
% integrate the kinematics in the extended window
opt = odeset('AbsTol',2.5e-14,'RelTol',2.3e-14);
[~,yy_ext] = ode113( @(t,y) kinematics_eom(t,y,J), tof_span_extended, [om0;q0],opt);
om_vec_ext = yy_ext(:,1:3)';
q_vec_ext  = yy_ext(:,4:7)';
% rotation matrices:
Rot_in2body = quat2dcm(q_vec');                % rotation matrix from inertial frame to body frame. [3x3xlength(t_span)] array.
Rot_in2lvlh = @(t) [ cos(n*t), sin(n*t), 0;
                    -sin(n*t), cos(n*t), 0;
                        0,         0,     1];  % rotation matrix from inertial to SGN LVLH frame.
% initial state estimate:
x0_mean = [ 15.792658268071492;       % initial position estimate [m] 
           -59.044939772661586; 
            3.227106250277039];
v0_mean = [-0.053960274403210;       % initial velocity estimate [m/s]
           -0.053969644762889;
           -0.089140748762173];
state0_mean = [x0_mean;v0_mean];
P0 = diag([10,10,10,0.1,0.1,0.1]);   % initial covariance matrix [m^2,m^2/s,m^2/s^2]
% nominal trajectory:
x0_nom = [12;-60;0];
v0_nom = [1e-4;-2*n*x0_nom(1);-1.2e-3];
state0_nom = [x0_nom;v0_nom];

%% camera model:
foc=30;  % [mm]
dens=54; % [pix/mm]
b=1;     % [m]
sens_size = [1920;1200];            % u0,v0 [pix]
p0=[sens_size(1)/2;sens_size(2)/2]; % center pixel location;
Cframe=[1,0,0;
        0,cos(-pi/2),sin(-pi/2);
        0,-sin(-pi/2),cos(-pi/2)];  % rotation from LVLH to camera frame.
R=10;
Cam.f=foc;
Cam.d=dens;
Cam.p0=p0;
Cam.b=b;
Cam.Cframe=Cframe;
Cam.R=R;
%% EX3.1
% state propagation:
fprintf('Propagating the states.\n');
tic;
% initial estimated mean state along 1 day window:
states_lvlh = CW_state_propagation(n,tof_span,state0_mean); 
% initial estimated mean state along 3T window:
states_lvlh_ext = CW_state_propagation(n,tof_span_extended,state0_mean);
% initial nominal state along 1 day window:
states_lvlh_nom = CW_state_propagation(n,tof_span,state0_nom);
% initial nominal state along 3T window:
states_lvlh_nom_ext = CW_state_propagation(n,tof_span_extended,state0_nom);
toc;
%% estimated state plots
figure();
hold on;
grid on;
plot3(0,0,0,'.k','LineWidth',3.5,'MarkerSize',30);
plot3(states_lvlh_ext(1,:)/1e3,states_lvlh_ext(2,:)/1e3,states_lvlh_ext(3,:)/1e3,'b','LineWidth',3);
plot3(states_lvlh(1,:)/1e3,states_lvlh(2,:)/1e3,states_lvlh(3,:)/1e3,'r','LineWidth',3.5);
xlabel('$x$ [km]'); ylabel('$y$[km]'), zlabel('$z$ [km]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15; ax.ZAxis.FontSize = 15;
legend('SGN-1 ',['$',num2str(n_periods),'$ orbital revolutions'],'One day navigation window','FontSize',15);
title('Relative orbit in SGN-1 LVLH frame using initial estimated state','FontSize',15);

figure();
hold on;
grid on;
plot(tof_span_extended/3600,states_lvlh_ext(1,:)/1e3,'b','LineWidth',3);
plot(tof_span/3600,states_lvlh(1,:)/1e3,'r','LineWidth',3);
xlabel('Time after $t_0$ [h]'); ylabel('$x$ [km]');
xlim([tof_span(1)/3600 tof_span_extended(end)/3600]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend(['$',num2str(n_periods),'$ orbital revolutions'],'One day navigation window','FontSize',15);
title('Relative orbit in SGN-1 LVLH frame using initial estimated state','FontSize',15);

figure();
hold on;
grid on;
plot(tof_span_extended/3600,states_lvlh_ext(2,:)/1e3,'b','LineWidth',3);
plot(tof_span/3600,states_lvlh(2,:)/1e3,'r','LineWidth',3);
xlabel('Time after $t_0$ [h]'); ylabel('$y$ [km]');
xlim([tof_span(1)/3600 tof_span_extended(end)/3600]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend(['$',num2str(n_periods),'$ orbital revolutions'],'One day navigation window','FontSize',15);
title('Relative orbit in SGN-1 LVLH frame using initial estimated state','FontSize',15);

figure();
hold on;
grid on;
plot(tof_span_extended/3600,states_lvlh_ext(3,:)/1e3,'b','LineWidth',3);
plot(tof_span/3600,states_lvlh(3,:)/1e3,'r','LineWidth',3);
xlabel('Time after $t_0$ [h]'); ylabel('$z$ [km]');
xlim([tof_span(1)/3600 tof_span_extended(end)/3600]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend(['$',num2str(n_periods),'$ orbital revolutions'],'One day navigation window','FontSize',15);
title('Relative orbit in SGN-1 LVLH frame using initial estimated state','FontSize',15);

%% nominal trajectory plots
figure();
hold on;
grid on;
plot3(states_lvlh_nom(1,:),states_lvlh_nom(2,:),states_lvlh_nom(3,:),'r','LineWidth',3);
plot3(0,0,0,'.b','LineWidth',3.5,'MarkerSize',30);
xlabel('$x$ [m]'); ylabel('$y$ [m]'), zlabel('$z$ [m]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15; ax.ZAxis.FontSize = 15;
legend('OSR nominal orbit','SGN-1','FontSize',15);
title('Nominal relative orbit in SGN-1 LVLH frame','FontSize',15);

figure();
hold on;
grid on;
plot(tof_span_extended/3600,states_lvlh_nom_ext(1,:),'b','LineWidth',3);
plot(tof_span/3600,states_lvlh_nom(1,:),'r','LineWidth',3);
xlabel('Time after $t_0$ [h]'); ylabel('$x$ [m]');
xlim([tof_span(1)/3600 tof_span_extended(end)/3600]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend(['$',num2str(n_periods),'$ orbital revolutions'],'One day navigation window','FontSize',15);
title('Nominal relative orbit in SGN-1 LVLH frame','FontSize',15);

figure();
hold on;
grid on;
plot(tof_span_extended/3600,states_lvlh_nom_ext(2,:),'b','LineWidth',3);
plot(tof_span/3600,states_lvlh_nom(2,:),'r','LineWidth',3);
xlabel('Time after $t_0$ [h]'); ylabel('$y$ [m]');
xlim([tof_span(1)/3600 tof_span_extended(end)/3600]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend(['$',num2str(n_periods),'$ orbital revolutions'],'One day navigation window','FontSize',15);
title('Nominal relative orbit in SGN-1 LVLH frame','FontSize',15);

figure();
hold on;
grid on;
plot(tof_span_extended/3600,states_lvlh_nom_ext(3,:),'b','LineWidth',3);
plot(tof_span/3600,states_lvlh_nom(3,:),'r','LineWidth',3);
xlabel('Time after $t_0$ [h]'); ylabel('$z$ [m]');
xlim([tof_span(1)/3600 tof_span_extended(end)/3600]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend(['$',num2str(n_periods),'$ orbital revolutions'],'One day navigation window','FontSize',15);
title('Nominal relative orbit in SGN-1 LVLH frame','FontSize',15);

figure();
hold on;
grid on;
plot(tof_span_extended/3600,states_lvlh_nom_ext(4,:),'b','LineWidth',3);
plot(tof_span/3600,states_lvlh_nom(4,:),'r','LineWidth',3);
xlabel('Time after $t_0$ [h]'); ylabel('$v_x \left[\frac{m}{s}\right]$');
xlim([tof_span(1)/3600 tof_span_extended(end)/3600]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend(['$',num2str(n_periods),'$ orbital revolutions'],'One day navigation window','FontSize',15);
title('Nominal relative orbit in SGN-1 LVLH frame','FontSize',15);

figure();
hold on;
grid on;
plot(tof_span_extended/3600,states_lvlh_nom_ext(5,:),'b','LineWidth',3);
plot(tof_span/3600,states_lvlh_nom(5,:),'r','LineWidth',3);
xlabel('Time after $t_0$ [h]'); ylabel('$v_y \left[\frac{m}{s}\right]$');
xlim([tof_span(1)/3600 tof_span_extended(end)/3600]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend(['$',num2str(n_periods),'$ orbital revolutions'],'One day navigation window','FontSize',15);
title('Nominal relative orbit in SGN-1 LVLH frame','FontSize',15);

figure();
hold on;
grid on;
plot(tof_span_extended/3600,states_lvlh_nom_ext(6,:),'b','LineWidth',3);
plot(tof_span/3600,states_lvlh_nom(6,:),'r','LineWidth',3);
xlabel('Time after $t_0$ [h]'); ylabel('$v_z \left[\frac{m}{s}\right]$');
xlim([tof_span(1)/3600 tof_span_extended(end)/3600]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend(['$',num2str(n_periods),'$ orbital revolutions'],'One day navigation window','FontSize',15);
title('Nominal relative orbit in SGN-1 LVLH frame','FontSize',15);

%% numbers of visible corners from the "real" sensor along nominal trajectory
fprintf('Computing real measurements in the nominal navigation window.\n');
tic;
[real_meas,n_vis_vertexes] = visibility_measurements(n,states_lvlh_nom(1:3,:),q_vec,tof_span,t0_ref,Cam);
toc;

figure();
scatter(tof_span/3600,n_vis_vertexes,25,'ob','filled');
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('Number of visible vertexes');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;

% compute trajectory of vertixes obtained from real measurements
[pos_vert1,tof_vert1] = real_vert_pos(1,real_meas,tof_span);
[pos_vert2,tof_vert2] = real_vert_pos(2,real_meas,tof_span);
[pos_vert3,tof_vert3] = real_vert_pos(3,real_meas,tof_span);
[pos_vert4,tof_vert4] = real_vert_pos(4,real_meas,tof_span);
[pos_vert5,tof_vert5] = real_vert_pos(5,real_meas,tof_span);
[pos_vert6,tof_vert6] = real_vert_pos(6,real_meas,tof_span);
[pos_vert7,tof_vert7] = real_vert_pos(7,real_meas,tof_span);
[pos_vert8,tof_vert8] = real_vert_pos(8,real_meas,tof_span);
%% EX3.2 vertex position using ideal model without accounting for visibility
N = length(tof_span);
sim_meas_vec = zeros(3,8,N);  % store results in a [3x8xN] matrix where third dimensions spans along the time epochs.
fprintf('Simulating measurements in nominal navigation window.\n');
tic;
for i = 1:N
    sim_meas_vec(:,:,i) = measurement_model(dim,states_lvlh_nom(1:3,i),q_vec(:,i),tof_span(i),Cam,Rot_in2lvlh);
end
toc;

% compare simulated and real measurements for each vertex:
figure();
hold on;
grid on;
plot(pos_vert1(1,:),pos_vert1(2,:),'.','Color','#EDB120','MarkerSize',15);
plot(squeeze(sim_meas_vec(1,1,:)),squeeze(sim_meas_vec(2,1,:)),'b','LineWidth',3);
xlim([-10 1400]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlabel('$q_{ul}$ [pix]','FontSize',15);ylabel('$q_{vl}$ [pix]','FontSize',15);
legend('Real measurements with visibility','Simulated measurements','FontSize',15);
title('Comparison of vertex $1$ trajectories','FontSize',15);

figure();
hold on;
grid on;
plot(pos_vert2(1,:),pos_vert2(2,:),'.','Color','#EDB120','MarkerSize',15);
plot(squeeze(sim_meas_vec(1,2,:)),squeeze(sim_meas_vec(2,2,:)),'b','LineWidth',3);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlabel('$q_{ul}$ [pix]','FontSize',15);ylabel('$q_{vl}$ [pix]','FontSize',15);
legend('Real measurements with visibility','Simulated measurements','FontSize',15);
title('Comparison of vertex $2$ trajectories','FontSize',15);

figure();
hold on;
grid on;
plot(pos_vert3(1,:),pos_vert3(2,:),'.','Color','#EDB120','MarkerSize',15);
plot(squeeze(sim_meas_vec(1,3,:)),squeeze(sim_meas_vec(2,3,:)),'b','LineWidth',3);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlabel('$q_{ul}$ [pix]','FontSize',15);ylabel('$q_{vl}$ [pix]','FontSize',15);
legend('Real measurements with visibility','Simulated measurements','FontSize',15);
title('Comparison of vertex $3$ trajectories','FontSize',15);

figure();
hold on;
grid on;
plot(pos_vert4(1,:),pos_vert4(2,:),'.','Color','#EDB120','MarkerSize',15);
plot(squeeze(sim_meas_vec(1,4,:)),squeeze(sim_meas_vec(2,4,:)),'b','LineWidth',3);
xlim([-10 1400]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlabel('$q_{ul}$ [pix]','FontSize',15);ylabel('$q_{vl}$ [pix]','FontSize',15);
legend('Real measurements with visibility','Simulated measurements','FontSize',15);
title('Comparison of vertex $4$ trajectories','FontSize',15);

figure();
hold on;
grid on;
plot(pos_vert5(1,:),pos_vert5(2,:),'.','Color','#EDB120','MarkerSize',15);
plot(squeeze(sim_meas_vec(1,5,:)),squeeze(sim_meas_vec(2,5,:)),'b','LineWidth',3);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlabel('$q_{ul}$ [pix]','FontSize',15);ylabel('$q_{vl}$ [pix]','FontSize',15);
legend('Real measurements with visibility','Simulated measurements','FontSize',15);
title('Comparison of vertex $5$ trajectories','FontSize',15);

figure();
hold on;
grid on;
plot(pos_vert6(1,:),pos_vert6(2,:),'.','Color','#EDB120','MarkerSize',15);
plot(squeeze(sim_meas_vec(1,6,:)),squeeze(sim_meas_vec(2,6,:)),'b','LineWidth',3);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlabel('$q_{ul}$ [pix]','FontSize',15);ylabel('$q_{vl}$ [pix]','FontSize',15);
legend('Real measurements with visibility','Simulated measurements','FontSize',15);
title('Comparison of vertex $6$ trajectories','FontSize',15);

figure();
hold on;
grid on;
plot(pos_vert7(1,:),pos_vert7(2,:),'.','Color','#EDB120','MarkerSize',15);
plot(squeeze(sim_meas_vec(1,7,:)),squeeze(sim_meas_vec(2,7,:)),'b','LineWidth',3);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlabel('$q_{ul}$ [pix]','FontSize',15);ylabel('$q_{vl}$ [pix]','FontSize',15);
legend('Real measurements with visibility','Simulated measurements','FontSize',15);
title('Comparison of vertex $7$ trajectories','FontSize',15);

figure();
hold on;
grid on;
plot(pos_vert8(1,:),pos_vert8(2,:),'.','Color','#EDB120','MarkerSize',15);
plot(squeeze(sim_meas_vec(1,8,:)),squeeze(sim_meas_vec(2,8,:)),'b','LineWidth',3);
xlim([-10 1400]);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
xlabel('$q_{ul}$ [pix]','FontSize',15);ylabel('$q_{vl}$ [pix]','FontSize',15);
legend('Real measurements with visibility','Simulated measurements','FontSize',15);
title('Comparison of vertex $8$ trajectories','FontSize',15);

%% EX3.3 EKF
fprintf('EKF during nominal navigation window.\n');
tic;
[ekf_state,ekf_P,n_correction] = EKF(state0_mean,P0,n,dim,real_meas,q_vec,tof_span,Cam,Rot_in2lvlh);
toc;

index_non_vis = n_vis_vertexes == 0;
index_vis = ~index_non_vis;
tof_non_vis = tof_span(index_non_vis);
tof_vis = tof_span(index_vis);
ekf_state_vis = ekf_state(:,index_vis);
ekf_state_non_vis = ekf_state(:,index_non_vis);

% plot ekf results in terms of error:

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(1,:) - ekf_state(1,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ekf_P(1,1,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$x$ [m]');
legend('$|x_{NOM} - x_{EKF}|$','$3\sqrt{P_{EKF}\left(1,1\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$x$ coordinate error using the EKF','FontSize',15);

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(2,:) - ekf_state(2,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ekf_P(2,2,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$y$ [m]');
legend('$|y_{NOM} - y_{EKF}|$','$3\sqrt{P_{EKF}\left(2,2\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$y$ coordinate error using the EKF','FontSize',15);

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(3,:) - ekf_state(3,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ekf_P(3,3,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$z$ [m]');
legend('$|z_{NOM} - z_{EKF}|$','$3\sqrt{P_{EKF}\left(3,3\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$z$ coordinate error using the EKF','FontSize',15);

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(4,:) - ekf_state(4,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ekf_P(4,4,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$v_x \left[\frac{m}{s}\right]$');
legend('$|v_{x_{NOM}} - v_{x_{EKF}}|$','$3\sqrt{P_{EKF}\left(4,4\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$v_x$ coordinate error using the EKF','FontSize',15);

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(5,:) - ekf_state(5,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ekf_P(5,5,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$v_y \left[\frac{m}{s}\right]$');
legend('$|v_{y_{NOM}} - v_{y_{EKF}}|$','$3\sqrt{P_{EKF}\left(5,5\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$v_y$ coordinate error using the EKF','FontSize',15);

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(6,:) - ekf_state(6,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ekf_P(6,6,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$v_z \left[\frac{m}{s}\right]$');
legend('$|v_{z_{NOM}} - v_{z_{EKF}}|$','$3\sqrt{P_{EKF}\left(6,6\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$v_z$ coordinate error using the EKF','FontSize',15);

%% EX3.3b UKF

% UKF in nominal navigation window
fprintf('UKF in nominal navigation window.\n');
tic;
[ukf_state,ukf_P,n_correction_ukf] = UKF(state0_mean,P0,n,dim,q_vec,real_meas,tof_span,Cam,Rot_in2lvlh);
toc;
ukf_state_vis = ukf_state(:,index_vis);
ukf_state_non_vis = ukf_state(:,index_non_vis);

% plot ukf results:
figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(1,:) - ukf_state(1,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ukf_P(1,1,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$x$ [m]');
legend('$|x_{NOM} - x_{UKF}|$','$3\sqrt{P_{UKF}\left(1,1\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$x$ coordinate error using the UKF','FontSize',15);

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(2,:)- ukf_state(2,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ukf_P(2,2,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$y$ [m]');
legend('$|y_{NOM} - y_{UKF}|$','$3\sqrt{P_{UKF}\left(2,2\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$y$ coordinate error using the UKF','FontSize',15);

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(3,:) - ukf_state(3,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ukf_P(3,3,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$z$ [m]');
legend('$|z_{NOM} - z_{UKF}|$','$3\sqrt{P_{UKF}\left(3,3\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$z$ coordinate error using the UKF','FontSize',15);

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(4,:) - ukf_state(4,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ukf_P(4,4,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$v_x \left[\frac{m}{s}\right]$');
legend('$|v_{x_{NOM}} - v_{x_{UKF}}|$','$3\sqrt{P_{UKF}\left(4,4\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$v_x$ coordinate error using the UKF','FontSize',15);

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(5,:) - ukf_state(5,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ukf_P(5,5,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$v_y \left[\frac{m}{s}\right]$');
legend('$|v_{y_{NOM}} - v_{y_{UKF}}|$','$3\sqrt{P_{UKF}\left(5,5\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$v_y$ coordinate error using the UKF','FontSize',15);

figure();
semilogy(tof_span/3600,abs(states_lvlh_nom(6,:) - ukf_state(6,:)),'b','LineWidth',2);
hold on;
semilogy(tof_span/3600,3*sqrt(squeeze(ukf_P(6,6,:))),'r','LineWidth',2);
grid on;
xlabel('Time after $t_0$ [h]'); ylabel('$v_z \left[\frac{m}{s}\right]$');
legend('$|v_{z_{NOM}} - v_{z_{UKF}}|$','$3\sqrt{P_{UKF}\left(6,6\right)}$','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('$v_z$ coordinate error using the UKF','FontSize',15);

%% Ex 3 functions

function dy = kinematics_eom(~,y,J)
% function that describes the evolution of the  attitude kinematics parameters.
% We assume Euler's euqation for the angular velocity behaviour and quaternions as kinematics
% parameters.
%
% INPUT:
% t = time. Not required because autonomous system.
% y = [7x1] vector of kinematics quantities defined as (omega;quaternion).
%     Angular velocity components are expressed in the spacecraft body frame.
%     Quaternions express orientation of this frame.
% J = [3x3] inertial matrix. In the body frame it is a diagonal matrix.
%
% OUTPUT:
% dy = [7x1] vector containing the derivatives of the attitude parameters.

om = y(1:3);
q = y(4:end);
dy = zeros(7,1);
skew_m = [0, om(3), -om(2); -om(3), 0, om(1); om(2), -om(1), 0];
% components of omega_dot:
dy(1:3) = J\(-cross(om,J*om));
% define matrix ruling quaternion evolution:
Qm = [ 0 , -om'; om, skew_m];
% components of q_dot:
dy(4:end) = 0.5*Qm*q;
end


function state = CW_state_equation(n,tof,state0)
% Function that computes the state at time t given initial condition.
% Based on the analytical solution for the Clohessy-Wiltshire model.
% INPUT:
% n      = mean motion of the target long its trajectory, [rev/sec].
% t      = time at which solution is to be computed, in elapsed seconds from t0 = 0, [sec].
% state0 = [6x1] initial conditions as (pos0,vel0), [km, km/sec].
%
% OUTPUT:
% state = [6x1] state at time t given as (pos,vel), [km, km/sec].
A = [   4-3*cos(n.*tof),       0,      0,       1/n*sin(n.*tof),      2/n*(1-cos(n.*tof)),               0;
     6*(sin(n.*tof)-n.*tof),   1,      0,     -2/n*(1-cos(n.*tof)),   1/n*(4*sin(n.*tof)-3*n.*tof),      0;
            0,                 0,  cos(n.*tof),        0,                     0,                    1/n*sin(n.*tof);
      3*n*sin(n.*tof),         0 ,     0,        cos(n.*tof),           2*sin(n.*tof),                   0;
      -6*n*(1-cos(n.*tof)),    0,      0,      -2*sin(n.*tof),          4*cos(n.*tof)-3,                 0;
           0,                  0, -n*sin(n.*tof),     0,                     0,                         cos(n.*tof)];
state = A*state0;
end

function states = CW_state_propagation(n,tof_span,state0)
% Function that computes the state at time t given initial condition.
% Based on the analytical solution for the Clohessy-Wiltshire model.
% INPUT:
% n        = mean motion of the target long its trajectory, [rev/sec].
% tof_span = [1xN] array containing times at which solution is to be 
% %                computed, defined as elapsed seconds from t0 = 0 [sec].
% state0   = [6x1] initial conditions as (pos0,vel0), [km, km/sec].
%
% OUTPUT:
% state = [6xN] state at times t_span given as (pos,vel), [km, km/sec].
N = length(tof_span);
states = zeros(length(state0),N);
for i = 1:N
    states(:,i) = CW_state_equation(n,tof_span(i),state0);
end
end


function [meas,n_vis_corners] = visibility_measurements(n,r_vec,q_vec,tof_vec,epoch0,Cam)
% function that extracts measurements seen by the camera at each time istant of t_vec.
% Relies on the usage of meas_sim_pvt.
% INPUT:
% n        = scalar representing the mean motion of the target spacecraft [rad/s]
% r        = [3xN] matrix representing the relative position of the chaser expressed in LVLH frame of the target @target [m].
% q_vec    = [4xN] matrix of quaternions with unit norm representing the attitude of the target w.r.t. inertial frame [-].
% tof_vec  = [1xN] vector representing the current instant of time passed from t0 [s]
% epoch0   = scalar representing the initial reference epoch expressed in TDB seconds past the J2000 (i.e., by using cspice_str2et) [s]
% Cam      = struct with the following fields:
%             - f       scalar that represents the focal length [mm]
%             - d       scalar that represents the pixel density [pix/mm]
%             - p0      2x1 vector containing the coordinates of the center pixel (u0;v0) [pix]
%             - b       scalar representing the stereo camera baseline [m]
%             - Cframe  3x3 director cosine matrix representing the rotation necessary to express a vector in LVLH frame to camera frame [-]
%             - R       scalar representing the variance of the measurement noise [pix]^2
%
% OUTPUT:
% meas          =  [1xN] struct array with the following fields:
%                 - y       [3xm] array of measurements containing the triplet of horizontal pixel, vertical pixel, and disparity for each of the m visible features
%                 - visible [1xm] array of IDs: each i-th entry identifies the vertex that generated the corresponding y(:,i) measurement  
% n_vis_corners =  [1xN] array containing numer of visible corners at each time istant
N = length(tof_vec);
meas = repmat(struct('y',0,'visible',0),1,N);
n_vis_corners = zeros(1,N);
for i = 1:N
    meas(i) = meas_sim_pvt(n,r_vec(:,i),q_vec(:,i),tof_vec(i),epoch0,Cam);
    n_vis_corners(i) = length(meas(i).visible); 
end
end


function meas = measurement_model(dim,r,q,tof,Cam,fun_in2lvlh)
% Function that simulates the pixel position of the s/c vertexes in the chaser's camera frame.
% Model used is an ideal stereo camera one. S/c assumed as a parallelepiped.
% 
% INPUT:
% dim         = [3x1] vector containing the s/c dimensions as (l,h,d). [m]
% r           = [3x1] chaser relative position in the target's lvlh frame. [m]
% q           = [4x1] current quaternion expressing orientation of target body frame w.r.t. inertial one.
% tof         = time passed from the reference epoch of the problem. [s]
% Cam         = struct with the following fields:
%                  - f       scalar that represents the focal length [mm]
%                  - d       scalar that represents the pixel density [pix/mm]
%                  - p0      2x1 vector containing the coordinates of the center pixel (u0;v0) [pix]
%                  - b       scalar representing the stereo camera baseline [m]
%                  - Cframe  3x3 director cosine matrix representing the rotation necessary to express a vector in LVLH frame to camera frame [-]
%                  - R       scalar representing the variance of the measurement noise [pix^2]
% fun_in2lvlh = [3x3] rotation matrix function of time that expresses rotation from inertial to lvlh
%               frame.
%
% OUTPUT:
% meas = [3x8] matrix where columns represents horizontal pixel,vertical pixel,disparity of each of
%        the 8 vertexes.

l = dim(1);
h = dim(2);
d = dim(3);
foc =Cam.f;
dens = Cam.d;
u0 = Cam.p0(1);
v0 = Cam.p0(2);
b = Cam.b;
rot_in2lv = fun_in2lvlh(tof); % rotation from inertial to lvlh frame.
rot_in2bd = quat2dcm(q');     % rotation from inertial to target body frame.
rot_lv2cam = Cam.Cframe;      % rotation from lvlh to camera frame.
vert_bod = [ l/2, -d/2, -h/2;
             l/2,  d/2, -h/2;
             l/2,  d/2,  h/2;
             l/2, -d/2,  h/2;
            -l/2, -d/2. -h/2;
            -l/2,  d/2, -h/2;
            -l/2,  d/2,  h/2;
            -l/2, -d/2,  h/2]';          % [3x8] matrix expressing position of vertexes in target's body frame.
vert_lv = rot_in2lv*rot_in2bd'*vert_bod; % transform vert matrix from body to inertial frame and then from inertial to lvlh.
vert_lv_osr = vert_lv - r;               % change the origin from target to chaser.
vert_cam_osr = rot_lv2cam*vert_lv_osr;   % rotate in chaser centered, camera fixed frame.

% extract corresponding stereo-camera measurements:
x = vert_cam_osr(1,:);
y = vert_cam_osr(2,:);
z = vert_cam_osr(3,:);
meas  = [u0 - dens*foc*y./z; v0 + dens*foc*x./z; b*dens*foc./z];  % [3x8] matrix containing camera measurements
end



function [vert_pos,vis_tof] = real_vert_pos(index,meas,tof_span)
% function that computes all position of the vertex -i observed using real measurements.
%
% INPUT:
% index    = scalar between 1 and 8 that specifies the vertex of interest.
% meas     =  [1xN] struct array with the following fields:
%              - y       [3xm] array of measurements containing the triplet of horizontal pixel, vertical pixel, and disparity for each of the m visible features
%              - visible [1xm] array of IDs: each i-th entry identifies the vertex that generated the corresponding y(:,i) measurement   
% tof_span = [1xN] vector of delta times at which meas is computed.
%
% OUTPUT:
% vert_pos = [3xM] matrix containing the visible positions of vertex -i in the camera frame.
% vis_tof  = [1xM] vector of delta times at which the vertex is visible.

% initialize quantities
N = length(tof_span);
vert_pos = NaN(3,N);
vis_tof = NaN(1,N);
% retrieve positions of vertex -i
for i = 1:N
    m = length(meas(i).visible);
    for j = 1:m
        if meas(i).visible(j) == index
            vert_pos(:,i) = meas(i).y(:,j);
            vis_tof(i) = tof_span(i);
        end
    end
end
% shrink array retrieving only visible times
i_visibility = isnan(vis_tof);
i_visibility = ~i_visibility;
vis_tof = vis_tof(i_visibility);
vert_pos = vert_pos(:,i_visibility);
end

function STM = CW_STM(n,t_i,t_f)
% function that computes the State Transition Matrix between initial and final times
% of the CW dynamics equations. Analytical expression is available.
% 
% INPUT:
% n   = scalar representing the mean motion of the target spacecraft [rad/s]
% t_i = initial time. [s]
% t_f = final time. [s]
%
% OUTPUT:
% STM = [6x6]state transition matrix between t_i and t_f.
tau = t_f - t_i;
STM = [   4-3*cos(n.*tau),        0,      0,       1/n*sin(n.*tau),      2/n*(1-cos(n.*tau)),               0;
        6*(sin(n.*tau)-n.*tau),   1,      0,     -2/n*(1-cos(n.*tau)),   1/n*(4*sin(n.*tau)-3*n.*tau),      0;
               0,                 0,  cos(n.*tau),        0,                     0,                    1/n*sin(n.*tau);
         3*n*sin(n.*tau),         0,      0,        cos(n.*tau),           2*sin(n.*tau),                   0;
       -6*n*(1-cos(n.*tau)),      0,      0,      -2*sin(n.*tau),         4*cos(n.*tau)-3,                   0;
               0,                 0, -n*sin(n.*tau),     0,                     0,                        cos(n.*tau)];
end

function Jac = meas_Jac(dim,r,q,tof,Cam,fun_in2lvlh)
% Function that computes the jacobian of the measurement model at given epoch with respect to
% the SGN position. This is achieved by using tha chain rule: dh/dx = dh/dX*dX/dx.
% 
% INPUT:
% dim         = [3x1] vector containing the s/c dimensions as (l,h,d). [m]
% r           = [3x1] chaser relative position in the target's lvlh frame. [m]
% q           = [4x1] current quaternion expressing orientation of target body frame w.r.t. inertial one.
% tof         = time passed from the reference epoch of the problem. [s]
% Cam         = struct with the following fields:
%                  - f       scalar that represents the focal length [mm]
%                  - d       scalar that represents the pixel density [pix/mm]
%                  - p0      2x1 vector containing the coordinates of the center pixel (u0;v0) [pix]
%                  - b       scalar representing the stereo camera baseline [m]
%                  - Cframe  3x3 director cosine matrix representing the rotation necessary to express a vector in LVLH frame to camera frame [-]
%                  - R       scalar representing the variance of the measurement noise [pix^2]
% fun_in2lvlh = [3x3] rotation matrix function of time that expresses rotation from inertial to lvlh
%               frame.
%
% OUTPUT:
% Jac = [24x6] jacobian matrixof the measurement equation. We have a 3x6 matrix for each vertex
% position.

l = dim(1);
h = dim(2);
d = dim(3);
foc =Cam.f;
dens = Cam.d;
b = Cam.b;
rot_in2lv = fun_in2lvlh(tof); % rotation from inertial to lvlh frame.
rot_in2bd = quat2dcm(q');     % rotation from inertial to target body frame.
rot_lv2cam = Cam.Cframe;      % rotation from lvlh to camera frame.
% Given the s/c state (position and velocity of target in its lvlh frame) 
% we need the corresponding position of each vertex.
vert_bod = [ l/2, -d/2, -h/2;
             l/2,  d/2, -h/2;
             l/2,  d/2,  h/2;
             l/2, -d/2,  h/2;
            -l/2, -d/2. -h/2;
            -l/2,  d/2, -h/2;
            -l/2,  d/2,  h/2;
            -l/2, -d/2,  h/2]';          % [3x8] matrix expressing position of vertexes in target's body frame.
vert_lv = rot_in2lv*rot_in2bd'*vert_bod; % transform vert matrix from body to inertial frame and then from inertial to lvlh.
vert_lv_osr = vert_lv - r;               % change the origin from target to chaser.
vert_cam_osr = rot_lv2cam*vert_lv_osr;   % rotate in chaser centered, camera fixed frame.

% build the 24x6 jacobian matrix
Jac_fun = @(x) [0,             -dens*foc/x(3),      dens*foc*x(2)/(x(3)^2);
                dens*foc/x(3),       0,            -dens*foc*x(1)/(x(3)^2);
                0,                   0,              -b*dens*foc/(x(3)^2)];
Jac = -[Jac_fun(vert_cam_osr(:,1))*rot_lv2cam, zeros(3);
        Jac_fun(vert_cam_osr(:,2))*rot_lv2cam, zeros(3);
        Jac_fun(vert_cam_osr(:,3))*rot_lv2cam, zeros(3);
        Jac_fun(vert_cam_osr(:,4))*rot_lv2cam, zeros(3);
        Jac_fun(vert_cam_osr(:,5))*rot_lv2cam, zeros(3);
        Jac_fun(vert_cam_osr(:,6))*rot_lv2cam, zeros(3);
        Jac_fun(vert_cam_osr(:,7))*rot_lv2cam, zeros(3);
        Jac_fun(vert_cam_osr(:,8))*rot_lv2cam, zeros(3)];
end

function [ekf_state,ekf_P,n_correction] = EKF(state0,P0,n,dim,real_meas,q_vec,tof_span,Cam,fun_in2lvlh)
% Extended Kalman Filter procedure to refine the OSR state estimation along its trajectory
% given real measurements and expected ones.
%
% INPUT:
% state0      = [6x1] array of estimated initial state.
% P0          = [6x6] initial covariance matrix of the state.
% n           = scalar representing the mean motion of the target spacecraft [rad/s]
% dim         = [3x1] vector containing the s/c dimensions as (l,h,d). [m]
% real_meas   = [1xN] struct array of real measurements containing following fields for each epoch:
%                 - y       [3xm] array of measurements containing the triplet of horizontal pixel, vertical pixel, and disparity for each of the m visible features
%                 - visible [1xm] array of IDs: each i-th entry identifies the vertex that generated the corresponding y(:,i) measurement 
% q_vec       = [4xN] matrix containing the SGN quaternions at each time epoch.
% tof_span    = [1xN] array of elapsed time since the reference epoch. [s]
% Cam         = struct with the following fields:
%                  - f       scalar that represents the focal length [mm]
%                  - d       scalar that represents the pixel density [pix/mm]
%                  - p0      2x1 vector containing the coordinates of the center pixel (u0;v0) [pix]
%                  - b       scalar representing the stereo camera baseline [m]
%                  - Cframe  3x3 director cosine matrix representing the rotation necessary to express a vector in LVLH frame to camera frame [-]
%                  - R       scalar representing the variance of the measurement noise [pix^2]
% fun_in2lvlh = [3x3] rotation matrix function of time that expresses rotation from inertial to lvlh
%               frame.
%
% OUTPUT:
% mean_state = [6x2N-1] matrix of mean state estimate. 2N-1 length since at each epoch we have both
%              the a priori and the updated information, except for the initial state.
% P          = [6x6x2N-1] covariance matrix where the third dimension spans the time epochs.

N = length(tof_span);
ekf_state = zeros(length(state0),N);
ekf_P = zeros(length(P0),length(P0),N);
n_correction = 0;

% try and update the first point we have:
STM = CW_STM(n,tof_span(1),tof_span(1));
P_km = STM*P0*STM';
if isempty(real_meas(1).visible)
    % can't update values because target not visible
    x_kp = state0;
    P_kp = P_km;
    ekf_state(:,1) = x_kp;
    ekf_P(:,:,1) = P_kp;
else
    % values can be updated
    y_km = measurement_model(dim,state0(1:3),q_vec(:,1),tof_span(1),Cam,fun_in2lvlh);
    y_km = y_km(:,real_meas(1).visible);
    y_km = y_km(:);
    real_y = real_meas(1).y(:);
    Hk = meas_Jac(dim,state0(1:3),q_vec(:,1),tof_span(:,1),Cam,fun_in2lvlh);
    % retrieve only visible elements of the jacobian
    index = real_meas(1).visible;
    jac_index = zeros(1,3*length(index));
    for j = 1:length(index)
        jac_index(3*j-2:3*j) = 3*index(j)-2:3*index(j);
    end
    Hk = Hk(jac_index,:);   
    A = Hk*P_km*Hk' + diag(Cam.R*ones(1,length(y_km)));
    Kk = P_km*Hk'/A;
    x_kp = state0 + Kk*(real_y - y_km);
    P_kp = (eye(length(P_km)) - Kk*Hk)*P_km;
    ekf_state(:,1) = x_kp;
    ekf_P(:,:,1) = P_kp;
    n_correction = n_correction+1;
end

% start the filter's loop:
for k = 2:N
    x_km1 = ekf_state(:,k-1);
    P_km1 = ekf_P(:,:,k-1);
    % first propagation step: compute state at t = 0:
    x_k0 = CW_state_equation(n,-tof_span(k-1),x_km1);
    % second step: propagate from t = 0 to current elapsed time
    x_km = CW_state_equation(n,tof_span(k),x_k0);
    STM  = CW_STM(n,tof_span(k-1),tof_span(k));
    P_km = STM*P_km1*STM';
    % access real measurements and check target visibility
    if isempty(real_meas(k).visible)
        % target not visible so can't upgrade the state
        x_kp = x_km;
        P_kp = P_km;
        ekf_state(:,k) = x_kp;
        ekf_P(:,:,k) = P_kp;
    else
        n_correction = n_correction+1;
        % target is visible: simulate measurements and compare to real ones
        y_km = measurement_model(dim,x_km(1:3),q_vec(:,k),tof_span(k),Cam,fun_in2lvlh);   % meas at x_k
        Hk   = meas_Jac(dim,x_km(1:3),q_vec(:,k),tof_span(k),Cam,fun_in2lvlh);            % Jacobian at x_k
        % retrieve only visible ones:
        index = real_meas(k).visible;
        y_km = y_km(:,index);
        jac_index = zeros(1,3*length(index));
        for j = 1:length(index)
            jac_index(3*j-2:3*j) = 3*index(j)-2:3*index(j);
        end
        Hk = Hk(jac_index,:);
        %define noise matrix
        R = diag(Cam.R*ones(1,3*length(index))); % diagonal 3Mx3M matrix
        % compute Kalman gain
        Kk = (P_km*Hk')/(Hk*P_km*Hk' + R);
        % reshape real and simulated measurements:
        y_km = y_km(:);              % column vector
        y_real = real_meas(k).y(:);    % column vector
        % update the estimates:
        x_kp = x_km + Kk*(y_real - y_km);
        P_kp = (eye(length(P_km)) - Kk*Hk)*P_km;
        ekf_state(:,k) = x_kp;
        ekf_P(:,:,k) = P_kp;
    end
end
end




function [ukf_state,ukf_P,n_correction] = UKF(state0_mean,P0,mean_motion,dim,q_vec,real_meas,tof_span,Cam,fun_in2lvlh)
% Unscented Kalman Filter procedure to refine the OSR state estimation along its trajectory
% given real measurements and expected ones.
%
% INPUT:
% state0      = [6x1] array of estimated initial state.
% P0          = [6x6] initial covariance matrix of the state.
% n           = scalar representing the mean motion of the target spacecraft [rad/s]
% dim         = [3x1] vector containing the s/c dimensions as (l,h,d). [m]
% real_meas   = [1xN] struct array of real measurements containing following fields for each epoch:
%                 - y       [3xm] array of measurements containing the triplet of horizontal pixel, vertical pixel, and disparity for each of the m visible features
%                 - visible [1xm] array of IDs: each i-th entry identifies the vertex that generated the corresponding y(:,i) measurement 
% q_vec       = [4xN] matrix containing the SGN quaternions at each time epoch.
% tof_span    = [1xN] array of elapsed time since the reference epoch. [s]
% Cam         = struct with the following fields:
%                  - f       scalar that represents the focal length [mm]
%                  - d       scalar that represents the pixel density [pix/mm]
%                  - p0      2x1 vector containing the coordinates of the center pixel (u0;v0) [pix]
%                  - b       scalar representing the stereo camera baseline [m]
%                  - Cframe  3x3 director cosine matrix representing the rotation necessary to express a vector in LVLH frame to camera frame [-]
%                  - R       scalar representing the variance of the measurement noise [pix^2]
% fun_in2lvlh = [3x3] rotation matrix function of time that expresses rotation from inertial to lvlh
%               frame.
%
% OUTPUT:
% mean_state = [6x2N-1] matrix of mean state estimate. 2N-1 length since at each epoch we have both
%              the a priori and the updated information, except for the initial state.
% P          = [6x6x2N-1] covariance matrix where the third dimension spans the time epochs.

% initialize parameters
n_correction = 0;
N = length(tof_span);
n = length(state0_mean);
ukf_state = zeros(n,N);         % mean state vector
ukf_P = zeros(n,n,N);           % covariance matrix array

% set parameters for unscented transform:
alpha = 1e-3;
beta = 2;
k = 0;
c = alpha^2*(n + k);
W0_m = 1 - n/(alpha^2*(n + k));
Wi_m = 1/(2*alpha^2*(n + k));
W0_c = 2 - alpha^2 + beta - n/(alpha^2*(n + k));
Wi_c = Wi_m;
Wm_vec = [W0_m, Wi_m*ones(1,2*n)];
Wc_vec = [W0_c, Wi_c*ones(1,2*n)];

% update initial state with available measurements:
if isempty(real_meas(1).visible)
    ukf_state(:,1) = state0_mean;
    ukf_P(:,:,1) = P0;
else
    sqrtP = sqrtm(c*P0);
    sigma_points = [state0_mean, state0_mean + sqrtP, state0_mean - sqrtP];
    % obtain their expected measurements 
    gamma_points = zeros(3*length(real_meas(1).visible),2*n + 1);
    for i = 1:2*n + 1
        exp_meas = measurement_model(dim,sigma_points(1:3,i),q_vec(:,1),tof_span(1),Cam,fun_in2lvlh);
        exp_meas = exp_meas(:,real_meas(1).visible);
        gamma_points(:,i) = exp_meas(:);
    end
    x_km = sum(Wm_vec.*sigma_points,2);
    y_km = sum(Wm_vec.*gamma_points,2);
    P_km = Wc_vec.*(sigma_points - x_km)*(sigma_points - x_km)';
    Pee_k = Wc_vec.*(gamma_points - y_km)*(gamma_points - y_km)' + diag(Cam.R*ones(1,length(y_km)));
    Pxy_k = Wc_vec.*(sigma_points - x_km)*(gamma_points - y_km)';
    Kk = Pxy_k/Pee_k;
    real_y = real_meas(1).y;
    real_y = real_y(:);
    x_kp = x_km + Kk*(real_y - y_km);
    P_kp = P_km - Kk*Pee_k*Kk';
    ukf_state(:,1) = x_kp;
    ukf_P(:,:,1) = P_kp;
    n_correction = n_correction + 1;
end

% start the for loop
for k = 2:N
    x_km1 = ukf_state(:,k-1);
    P_km1 = ukf_P(:,:,k-1);
    % generate sigma points and propagate them:
    sqrtP = sqrtm(c*P_km1);
    sigma_points_km1 = [x_km1, x_km1 + sqrtP, x_km1 - sqrtP];
    sigma_points_k = zeros(n, 2*n + 1);
    % propagate sigma points
    for i = 1:2*n + 1
        % first propagate them backwards to -tof(k-1)
        sigma_points_k0 = CW_state_equation(mean_motion,-tof_span(k-1),sigma_points_km1(:,i));
        % then propagate forwards to tof(k)
        sigma_points_k(:,i) = CW_state_equation(mean_motion,tof_span(k),sigma_points_k0);
    end
    % compute mean and covariance
    x_km = sum(Wm_vec.*sigma_points_k,2);
    P_km = Wc_vec.*(sigma_points_k - x_km)*(sigma_points_k - x_km)';
    % update step
    if isempty(real_meas(k).visible)
        % update not possible
        x_kp = x_km;
        P_kp = P_km;
        ukf_state(:,k) = x_kp;
        ukf_P(:,:,k) = P_kp;
    else
        % update possible: obtain their expected measurements 
        gamma_points = zeros(3*length(real_meas(k).visible),2*n + 1);
        for i = 1:2*n + 1
            exp_meas = measurement_model(dim,sigma_points_k(1:3,i),q_vec(:,k),tof_span(k),Cam,fun_in2lvlh);
            exp_meas = exp_meas(:,real_meas(k).visible);
            gamma_points(:,i) = exp_meas(:);
        end
        % compute mean and covariance of predicted measurements:
        y_km = sum(Wm_vec.*gamma_points,2);
        Pee_k = Wc_vec.*(gamma_points - y_km)*(gamma_points - y_km)' + diag(Cam.R*ones(1,length(y_km)));
        Pxy_k = Wc_vec.*(sigma_points_k - x_km)*(gamma_points - y_km)';
        real_y = real_meas(k).y;
        real_y = real_y(:);
        % update step
        Kk = Pxy_k/Pee_k;
        x_kp = x_km + Kk*(real_y - y_km);
        P_kp = P_km - Kk*Pee_k*Kk';
        ukf_state(:,k) = x_kp;
        ukf_P(:,:,k) = P_kp;
        n_correction = n_correction + 1;
    end
end
end

