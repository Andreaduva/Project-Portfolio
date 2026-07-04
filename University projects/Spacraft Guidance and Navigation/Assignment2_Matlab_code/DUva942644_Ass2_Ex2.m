%% ASSIGNMENT 2 EX 2

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
%cspice_furnsh('..\assignment02_win.tm');  % Windows version
cspice_furnsh('../assignment02_mac.tm');  % MacOs version
kernel_number = cspice_ktotal('ALL');
fprintf('\nNumber of kernels loaded is %d\n',kernel_number);

%% INTRODUCTION
% set constants and parameters:
mu = cspice_bodvrd('EARTH','GM',1);
r_Earth = cspice_bodvrd('Earth','RADII',3);
t_ref0 = cspice_str2et('2022-11-11 19:08:49.824 UTC');
x0_ref = [6054.30795817484; -3072.03883303992; -133.115352431876];
v0_ref = [4.64750094824087; 9.18608475681236; -0.62056520749034];
state_ref0 = [x0_ref;v0_ref];
t0 = cspice_str2et('2022-11-12 04:30:00.000 UTC');
state_t0 = twobody_flow(state_ref0, 0, t0 - t_ref0, mu);
fprintf('Initial epoch t0 is 2022-11-12 04:30:00.000 UTC\n');
t_f = cspice_str2et('2022-11-14 16:30:00.000 UTC');
fprintf('Final epoch tf is 2022-11-14 16:30:00.000 UTC\n');
tof = t_f - t0;                                             % time of flight in seconds
tof_span = 0:60:tof;                                        % prediction every minute
et_vec = t0 + tof_span;
Kourou.name = 'KOUROU';
Perth.name  = 'PERTH';
Kourou.frame = 'KOUROU_TOPO';
Perth.frame = 'PERTH_TOPO';
Kourou.El_treshold = 10;          % deg
Perth.El_treshold = 5;            % deg

% plot Kourou and Perth reference frames:
% origin of Kourou and Perth:
Kourou.pos = cspice_spkpos('KOUROU',t0,'J2000','NONE','EARTH');
Perth.pos = cspice_spkpos('Perth',t0,'J2000','NONE','EARTH');
Eci2Kou = cspice_pxform('KOUROU_TOPO','J2000',t0);
Eci2Per = cspice_pxform('PERTH_TOPO','J2000',t0);
Kou_ijk = Eci2Kou*eye(3);
Per_ijk = Eci2Per*eye(3);
figure();
hold on;
grid on;
[X,Y,Z] = ellipsoid(0,0,0,r_Earth(1), r_Earth(2), r_Earth(3),50);
surf(X,Y,Z,'FaceColor','none');
quiver3(0,0,0,r_Earth(1)+1000,0,0,'LineWidth',4,'AutoScale','off');
quiver3(0,0,0,0,r_Earth(2)+1000,0,'LineWidth',4,'AutoScale','off');
quiver3(0,0,0,0,0,r_Earth(3)+1000,'LineWidth',4,'AutoScale','off');
plot3(Kourou.pos(1),Kourou.pos(2),Kourou.pos(3),'or','MarkerFaceColor','r','MarkerSize',10)
quiver3(Kourou.pos(1),Kourou.pos(2),Kourou.pos(3),Kou_ijk(1,1)'*4e3,Kou_ijk(2,1)'*4e3,Kou_ijk(3,1)'*4e3,'LineWidth',4);
quiver3(Kourou.pos(1),Kourou.pos(2),Kourou.pos(3),Kou_ijk(1,2)'*4e3,Kou_ijk(2,2)'*4e3,Kou_ijk(3,2)'*4e3,'LineWidth',4);
quiver3(Kourou.pos(1),Kourou.pos(2),Kourou.pos(3),Kou_ijk(1,3)'*4e3,Kou_ijk(2,3)'*4e3,Kou_ijk(3,3)'*4e3,'LineWidth',4);
xlabel('$x$ [km]'); ylabel('$y$ [km]'); zlabel('$z$ [km]');
legend('','$x$ axis ECI','$y$ axis ECI','$z$ axis ECI','Kourou station',...
    '$x$ axis Kourou-Topo','$y$ axis Kourou-Topo','$z$ axis Kourou-Topo','FontSize',15);
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15; ax.ZAxis.FontSize = 15;
axis equal;

%% EX2.1
% propagate the reference trajectory along the time interval from t_ref1 to t_ref2
opt_ode = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
[~,yy] = ode113(@(t,y) eom_2body(t,y,mu), tof_span, state_t0, opt_ode);

% obtain azimuth and elevation for both stations at each propagated epoch:
[~,Kourou.Az,Kourou.El] = ECI2RAzEl(Kourou.name,et_vec,yy');
[~,Perth.Az,Perth.El] = ECI2RAzEl(Perth.name,et_vec,yy');
% visibility windows:
[Kourou.et_vis,Kourou.Az_vis,Kourou.El_vis,Kourou.index] = visibility_window( [Kourou.Az; Kourou.El], et_vec, Kourou.El_treshold);
[Perth.et_vis,Perth.Az_vis,Perth.El_vis,Perth.index] = visibility_window( [Perth.Az; Perth.El], et_vec, Perth.El_treshold);

Kourou.index = double(Kourou.index);
Kourou.wind_index = 1:length(et_vec);
Kourou.wind_index = Kourou.wind_index(ischange(Kourou.index));
if Kourou.El(1) > Kourou.El_treshold               % first window already open so not detected
    Kourou.wind_index = [1,Kourou.wind_index];
end

Perth.index = double(Perth.index);
Perth.wind_index = 1:length(et_vec);
Perth.wind_index = Perth.wind_index(ischange(Perth.index));
if Perth.El(1) > Perth.El_treshold               % first window already open so not detected
    Perth.wind_index = [1,Perth.wind_index];
end

% find and print data for visibility windows
Kourou.vis_vec_index = ones(1,length(Kourou.wind_index));
for i = 1:length(Kourou.wind_index)         % loop on the indexes of the windows
    if mod(i,2) == 0                                 % if closing epoch
        Kourou.wind_index(i) = Kourou.wind_index(i)-1; % consider the previous one, current is out of the f.o.v.
        str = cspice_et2utc(et_vec(Kourou.wind_index(i)),'ISOC',4); % calendar format of epoch
        fprintf('End of Kourou visibility window %d occurs at epoch %s\n',i/2,str);
        Kourou.vis_vec_index(i) = Kourou.vis_vec_index(i-1) + Kourou.wind_index(i) - Kourou.wind_index(i-1);
    else                                             % opening epoch
        str = cspice_et2utc(et_vec(Kourou.wind_index(i)),'ISOC',4); % calendar format of epoch
        fprintf('Beginning of Kourou visibility window %d occurs at epoch %s\n',(i+1)/2,str);
        if i==1
        else
            Kourou.vis_vec_index(i) = Kourou.vis_vec_index(i-1) + 1;
        end
    end
end

Perth.vis_vec_index = ones(1,length(Perth.wind_index));
for i = 1:length(Perth.wind_index)         % loop on the indexes of the windows
    if mod(i,2) == 0                                 % if closing epoch
        Perth.wind_index(i) = Perth.wind_index(i)-1; % consider the previous one, current is out of the f.o.v.
        str = cspice_et2utc(et_vec(Perth.wind_index(i)),'ISOC',4); % calendar format of epoch
        fprintf('End of Perth visibility window %d occurs at epoch %s\n',i/2,str);
        Perth.vis_vec_index(i) = Perth.vis_vec_index(i-1) + Perth.wind_index(i) - Perth.wind_index(i-1);
    else                                             % opening epoch
        str = cspice_et2utc(et_vec(Perth.wind_index(i)),'ISOC',4); % calendar format of epoch
        fprintf('Beginning of Perth visibility window %d occurs at epoch %s\n',(i+1)/2,str);
        if i==1
        else
            Perth.vis_vec_index(i) = Perth.vis_vec_index(i-1) + 1;
        end
    end
end

% epochs associated to these windows:
Kourou.wind_et = et_vec(Kourou.wind_index);
Perth.wind_et = et_vec(Perth.wind_index);

% plot the visibility windows:
figure();
hold on; grid on;
plot(tof_span/3600,double(Kourou.index),'LineWidth',2);
plot(tof_span/3600,double(Perth.index),'LineWidth',2);
legend('Kourou','Perth','FontSize',15);
ylim([-0.5 1.5]);
xlabel('$t$ - $t_0$ [h]');
yticks([-0.5 0 1 1.5]);
yticklabels({'','Invisible','Visible',''});
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Visibility windows','FontSize',15);



% Elevation plots
figure();
grid on;
hold on;
plot(tof_span/3600,Kourou.El,'LineWidth',4);
yline(Kourou.El_treshold,'k','LineWidth',2);
ylim([-90 90]);
xlabel('$t$ - $t_0$ [h]');
ylabel('Elevation [Deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('Arian 5','Elevation treshold','FontSize',15);
title('Elevation @Kourou','FontSize',15);

figure();
grid on;
hold on;
plot(tof_span/3600,Perth.El,'Color',"#D95319",'LineWidth',4);
yline(Perth.El_treshold,'k','LineWidth',2);
ylim([-90 90]);
xlabel('$t$ - $t_0$ [h]');
ylabel('Elevation [deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('Arian 5','Elevation treshold','FontSize',15);
title('Elevation @Perth','FontSize',15);


% polar plots
figure();
hold on;
grid on;
for i = 1:length(Kourou.wind_index)
     if mod(i,2) == 0  % window closes
         label = ['Window number ',num2str(fix(i/2))];
         plot(Kourou.Az(Kourou.wind_index(i-1):Kourou.wind_index(i)),Kourou.El(Kourou.wind_index(i-1):Kourou.wind_index(i)),"*",'MarkerSize',3.5 ...
             ,'DisplayName',label);
     end
end
legend show; legend('FontSize',15);
xlim([-180 180]); ylim([Kourou.El_treshold 100]); 
xlabel('Azimuth [deg]');
ylabel('Elevation [deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Polar plot @Kourou','FontSize',15);


figure();
hold on;
grid on;
for i = 1:length(Perth.wind_index)
     if mod(i,2) == 0  % window closes
         label = ['Window number ',num2str(fix(i/2))];
         plot(Perth.Az(Perth.wind_index(i-1):Perth.wind_index(i)),Perth.El(Perth.wind_index(i-1):Perth.wind_index(i)),"*",'MarkerSize',3.5 ...
             ,'DisplayName',label);
     end
end
legend show; legend('FontSize',15);
xlim([-100 100]); ylim([Perth.El_treshold 60]); 
xlabel('Azimuth [deg]');
ylabel('Elevation [deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Polar plot @Perth','FontSize',15);

%% EX2.2a

% define parameters required for use of sgp4
arcsec2rad = pi / (180*3600);
typerun    = 'u';                    % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';                    % afspc approach ('air force space command')
whichconst =  72;                    % WGS72 constants (radius, gravitational parameter)
% call the read_3LE function
[sat_rec,string1,string2] = read_3LE('87654','orbits.3le',whichconst);

% Get TLE epoch
[year,mon,day,hr,min,sec] = invjday(sat_rec.jdsatepoch, sat_rec.jdsatepochf);
sgp4_epoch_str = sprintf('%d-%02d-%02d %02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);

sgp4_epoch_tle = cspice_str2et(sgp4_epoch_str);
fprintf('Satellite num ID: %d\n', sat_rec.satnum);
fprintf('TLE reference epoch: %s UTC\n', sgp4_epoch_str);

% correction values for nutation and precession at TLE reference epoch
dPsi = -0.112967*arcsec2rad;
dEps = -0.006958*arcsec2rad;

% propagate state using sgp4 for both stations in the visibility window
fprintf('Point 2.2a: begin orbit propagation using SGP4.\n');
tic;
[sgp4_rv_eci_Kvis(1:3,:), sgp4_rv_eci_Kvis(4:6,:)] = sgp4_propagator(sat_rec, Kourou.et_vis, dPsi, dEps);
[sgp4_rv_eci_Pvis(1:3,:), sgp4_rv_eci_Pvis(4:6,:)] = sgp4_propagator(sat_rec, Perth.et_vis, dPsi, dEps);
fprintf('Process took %f seconds.\n',toc);
% retrieve measured range, azimuth, elevation.
[Kourou.sgp4.range_vis, Kourou.sgp4.Az_vis, Kourou.sgp4.El_vis] = ECI2RAzEl(Kourou.name, Kourou.et_vis, sgp4_rv_eci_Kvis);
[Perth.sgp4.range_vis, Perth.sgp4.Az_vis, Perth.sgp4.El_vis] = ECI2RAzEl(Perth.name, Perth.et_vis, sgp4_rv_eci_Pvis);

figure();
grid on;
hold on;
plot((Kourou.et_vis-t0)/3600,abs(Kourou.sgp4.El_vis)- abs(Kourou.El_vis),'o','Color',"#0072BD",'MarkerSize',4,'MarkerFaceColor',"#0072BD");
xlabel('$t$ - $t_0$ [h]'); ylabel('$\Delta$El [Deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Difference in elevation @Kourou','FontSize',15);

figure();
grid on;
hold on;
plot((Kourou.et_vis-t0)/3600,abs(Kourou.sgp4.Az_vis) - abs(Kourou.Az_vis),'o','Color',"#D95319",'MarkerSize',4,'MarkerFaceColor',"#D95319");
xlabel('$t$ - $t_0$ [h]'); ylabel('$\Delta$Az [Deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Difference in Azimuth @Kourou','FontSize',15);

figure();
grid on;
hold on;
plot((Perth.et_vis-t0)/3600,abs(Perth.sgp4.El_vis) - abs(Perth.El_vis),'o','Color',"#0072BD",'MarkerSize',4,'MarkerFaceColor',"#0072BD");
xlabel('$t$ - $t_0$ [h]'); ylabel('$\Delta$El [Deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Difference in elevation @Perth','FontSize',15);

figure();
grid on;
hold on;
plot((Perth.et_vis-t0)/3600,abs(Perth.sgp4.Az_vis)-abs(Perth.Az_vis),'o','Color',"#D95319",'MarkerSize',4,'MarkerFaceColor',"#D95319");
xlabel('$t$ - $t_0$ [h]'); ylabel('$\Delta$Az [Deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Difference in azimuth @Perth','FontSize',15);

%% EX2.2b

% define noise matrix R which is equal for both stations
sigma_r = 0.01;   % km
sigma_AzEl = 0.1;  % deg
R = diag([sigma_r sigma_AzEl sigma_AzEl])^2;           % units are [km^2, deg^2, deg^2]
W_meas = diag(1./[sigma_r, sigma_AzEl, sigma_AzEl]);   % Actually equals to R*-1/2

Kourou.sgp4.measures = [Kourou.sgp4.range_vis',Kourou.sgp4.Az_vis', Kourou.sgp4.El_vis'];
Kourou.sgp4.measures_noise =  mvnrnd(Kourou.sgp4.measures,R);

Perth.sgp4.measures = [Perth.sgp4.range_vis',Perth.sgp4.Az_vis', Perth.sgp4.El_vis'];
Perth.sgp4.measures_noise =  mvnrnd(Perth.sgp4.measures,R);

% plot random behaviour of noise
figure();
hold on;
grid on;
for i = 1:length(Kourou.wind_index)
     if mod(i,2) == 0  % window closes
         label = ['Window ',num2str(fix(i/2))];
         plot( ( Kourou.et_vis(Kourou.vis_vec_index(i-1):Kourou.vis_vec_index(i)) - t0) /3600,...
             abs( Kourou.sgp4.measures( Kourou.vis_vec_index(i-1):Kourou.vis_vec_index(i),1 ) )...
             - abs( Kourou.sgp4.measures_noise( Kourou.vis_vec_index(i-1):Kourou.vis_vec_index(i),1 ) ), ...
             'DisplayName',label,'LineWidth',2);
     end
end
yline(3*sigma_r,'k','LineWidth',2,'HandleVisibility','off'); 
yline(-3*sigma_r,'k','LineWidth',2,'DisplayName','$\pm 3 \sigma_r$');
legend show; legend('FontSize',15); ylim([-0.08 0.08]);
xlabel('$t$ - $t_0$ [h]');
ylabel('range error [km]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Measurement noise @Kourou','FontSize',15);


figure();
hold on;
grid on;
for i = 1:length(Kourou.wind_index)
     if mod(i,2) == 0  % window closes
         label = ['Window ',num2str(fix(i/2))];
         plot( ( Kourou.et_vis(Kourou.vis_vec_index(i-1):Kourou.vis_vec_index(i)) - t0) /3600,...
             abs( Kourou.sgp4.measures( Kourou.vis_vec_index(i-1):Kourou.vis_vec_index(i),2 ) )...
             - abs( Kourou.sgp4.measures_noise( Kourou.vis_vec_index(i-1):Kourou.vis_vec_index(i),2 ) ), ...
             'DisplayName',label,'LineWidth',2);
     end
end
yline(3*sigma_AzEl,'k','LineWidth',2,'HandleVisibility','off'); 
yline(-3*sigma_AzEl,'k','LineWidth',2,'DisplayName','$\pm 3 \sigma_{Az}$');
legend show; legend('FontSize',15);ylim([-0.8 0.8]);
xlabel('$t$ - $t_0$ [h]');
ylabel('Azimuth error [deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Measurement noise @Kourou','FontSize',15);


figure();
hold on;
grid on;
for i = 1:length(Perth.wind_index)
     if mod(i,2) == 0  % window closes
         label = ['Window ',num2str(fix(i/2))];
         plot( ( Perth.et_vis(Perth.vis_vec_index(i-1):Perth.vis_vec_index(i)) - t0) /3600,...
             abs( Perth.sgp4.measures( Perth.vis_vec_index(i-1):Perth.vis_vec_index(i),3 ) )...
             - abs( Perth.sgp4.measures_noise( Perth.vis_vec_index(i-1):Perth.vis_vec_index(i),3 ) ), ...
             'DisplayName',label,'LineWidth',2);
     end
end
yline(3*sigma_AzEl,'k','LineWidth',2,'HandleVisibility','off'); 
yline(-3*sigma_AzEl,'k','LineWidth',2,'DisplayName','$\pm 3 \sigma_r$');
legend show; legend('FontSize',15);ylim([-0.8 0.8]);
xlabel('$t$ - $t_0$ [h]');
ylabel('Elevation error [deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Measurement noise @Perth','FontSize',15);

%% EX 2.3a
% reference epoch and state:
[pos0, vel0] = sgp4_propagator(sat_rec,t0,dPsi,dEps);
state0 = [pos0;vel0];
% transpose array for simplicity
Kourou.sgp4.measures = Kourou.sgp4.measures';
Kourou.sgp4.measures_noise = Kourou.sgp4.measures_noise';
Perth.sgp4.measures = Perth.sgp4.measures';
Perth.sgp4.measures_noise = Perth.sgp4.measures_noise';

% parameters for function that computes residuals
parameters_a.t_ref = t0;
parameters_a.et_vis = Perth.et_vis;
parameters_a.station_ID = 2*ones(length(Perth.et_vis),1);
parameters_a.mu = mu;
parameters_a.measurements = Perth.sgp4.measures_noise;
parameters_a.motion = 'Keplerian';
parameters_a.W_meas = W_meas;

% propagate refence state;
[~,yy_ref] = ode113(@(t,y) eom_2body(t,y,mu), et_vec - t0, state0, opt_ode);
yy_ref_a = yy_ref(logical(Perth.index),:)';
% compute associated measures:
[range_ref_a, Az_ref_a, El_ref_a] = ECI2RAzEl(Perth.name, Perth.et_vis, yy_ref_a);

% initialize residuals vector
fprintf('Point 2.3a: computing residuals with respect to reference trajectory.\n');
tic;
residual_ref_a = zeros(3,length(Perth.et_vis));   
% start a loop to compute residuals at each epoch:
for i = 1:length(Perth.et_vis)
    meas_pred = [range_ref_a(i); Az_ref_a(i); El_ref_a(i)];
    meas_real = Perth.sgp4.measures_noise(:,i);
    residual_ref_a(:,i) = W_meas*(meas_pred - meas_real);    % residual = weight*(pred_meas - real_meas)
end
fprintf('Process took %f seconds\n',toc);


% Solve nonlinear least squares problem
opt_lsq = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter-detailed');
fun_a = @(x) costfunction_radar(x,parameters_a);
%call solver
fprintf('Point 2.3a: optimization process using measurements at Perth and Keplerian motion.\n'); 
tic;
[x_sol_a, resnorm_a, residual_a, flag_a, ~, ~, Jac_a] = lsqnonlin(fun_a,state0, [], [], opt_lsq);
fprintf('Optimization took %f seconds.\n',toc);
% estimation of covariance matrix:
Jac_a = full(Jac_a);
P_ls_a = resnorm_a/( length(residual_a) - length(state0) ) .*inv(Jac_a'*Jac_a);

% propagate solution state;
[~,yy_sol] = ode113(@(t,y) eom_2body(t,y,mu), et_vec - t0, x_sol_a, opt_ode);
yy_sol_a = yy_sol(logical(Perth.index),:)';
% compute associated measures:
[range_sol_a, Az_sol_a, El_sol_a] = ECI2RAzEl(Perth.name, Perth.et_vis, yy_sol_a);
% compute residuals associated to the solution
residual_sol_a = zeros(3,length(Perth.et_vis));   
% start a loop to compute residuals at each epoch:
for i = 1:length(Perth.et_vis)
    meas_pred = [range_sol_a(i); Az_sol_a(i); El_sol_a(i)];
    meas_real = Perth.sgp4.measures_noise(:,i);
    residual_sol_a(:,i) = W_meas*(meas_pred - meas_real);    % residual = weight*(pred_meas - real_meas)
end
fprintf('Process took %f seconds\n',toc);

figure();
hold on;
grid on;
plot((Perth.et_vis-t0)/3600, residual_ref_a(1,:)./W_meas(1,1),'.b','MarkerSize',15,'MarkerFaceColor','b');
plot((Perth.et_vis-t0)/3600, residual_sol_a(1,:)./W_meas(1,1),'.r','MarkerSize',15,'MarkerFaceColor','b');
xlabel('$t$ - $t_0$ [h]'); ylabel('range [km]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('Reference residuals','Best estimate residuals','FontSize',15);
title('Perth measurements only and Keplerian motion','FontSize',15);

figure();
hold on;
grid on;
plot((Perth.et_vis-t0)/3600, residual_ref_a(2,:)./W_meas(2,2),'.b','MarkerSize',15,'MarkerFaceColor','b');
plot((Perth.et_vis-t0)/3600, residual_sol_a(2,:)./W_meas(2,2),'.r','MarkerSize',15,'MarkerFaceColor','b');
xlabel('$t$ - $t_0$ [h]'); ylabel('$Az$ [deg]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
legend('Reference residuals','Best estimate residuals','FontSize',15);
title('Perth measurements and Keplerian motion','FontSize',15);
%% EX 2.3b
% initialize global quantities:
et_vis_global = [Perth.et_vis,Kourou.et_vis];
meas_real_global = [Perth.sgp4.measures_noise, Kourou.sgp4.measures_noise];
station_ID_global = [2*ones(1,length(Perth.et_vis)), ones(1,length(Kourou.et_vis))];

% sort them by date:
[et_vis_global, sort_index] = sort(et_vis_global);
meas_real_global = meas_real_global(:,sort_index);
station_ID_global = station_ID_global(sort_index);

parameters_b.t_ref = t0;
parameters_b.et_vis = et_vis_global;
parameters_b.station_ID = station_ID_global;
parameters_b.mu = mu;
parameters_b.measurements = meas_real_global;
parameters_b.motion = 'Keplerian';
parameters_b.W_meas = W_meas;

% propagate reference state
if et_vis_global(1) - t0 == 0
    [~,yy_ref] = ode113(@(t,y) eom_2body(t,y,mu), et_vis_global - t0, state0, opt_ode);
    yy_ref_b = yy_ref(:,:)';
elseif et_vis_global(1) - t0 ~= 0
    [~,yy_ref] = ode113(@(t,y) eom_2body(t,y,mu), [0, et_vis_global - t0], state0, opt_ode);
    yy_ref_b = yy_ref(2:end,:)';
end

% initialize residuals vector
fprintf('Point 2.3b: computing residuals with respect to reference trajectory.\n');
tic;
residual_ref_b = zeros(3,length(et_vis_global));   
% start a loop to compute residuals at each epoch:
for i = 1:length(et_vis_global)
    if station_ID_global(i) == 1
        station_name = 'KOUROU';
    elseif station_ID_global(i) == 2
        station_name = 'PERTH';
    end
    [range_pred,Az_pred,El_pred] = ECI2RAzEl(station_name,et_vis_global(i),yy_ref_b(:,i));
    meas_pred = [range_pred; Az_pred; El_pred];
    meas_real = meas_real_global(:,i);
    residual_ref_b(:,i) = W_meas*(meas_pred - meas_real);    % residual = weight*(pred_meas - real_meas)
end
fprintf('Process took %f seconds\n',toc);


% Solve nonlinear least squares problem
opt_lsq = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter-detailed');
fun_b = @(x) costfunction_radar(x,parameters_b);
%call solver
fprintf('Point 2.3b: optimization process assuming Keplerian motion while using measurements at Perth and Kourou.\n'); 
tic;
[x_sol_b, resnorm_b, residual_b, flag_b, ~, ~, Jac_b] = lsqnonlin(fun_b,state0, [], [], opt_lsq);
fprintf('Optimization took %f seconds.\n',toc);
% estimation of covariance matrix:
Jac_b = full(Jac_b);
P_ls_b = resnorm_b/( length(residual_b) - length(state0) ) .*inv(Jac_b'*Jac_b);

% propagate solution state
if et_vis_global(1) - t0 == 0
    [~,yy_sol] = ode113(@(t,y) eom_2body(t,y,mu), et_vis_global - t0, x_sol_b, opt_ode);
    yy_sol_b = yy_sol(:,:)';
elseif et_vis_global(1) - t0 ~= 0
    [~,yy_sol] = ode113(@(t,y) eom_2body(t,y,mu), [0, et_vis_global - t0], x_sol_b, opt_ode);
    yy_sol_b = yy_sol(2:end,:)';
end

% initialize residuals vector
fprintf('Point 2.3b: computing residuals with respect to reference trajectory.\n');
tic;
residual_sol_b = zeros(3,length(et_vis_global));   
% start a loop to compute residuals at each epoch:
for i = 1:length(et_vis_global)
    if station_ID_global(i) == 1
        station_name = 'KOUROU';
    elseif station_ID_global(i) == 2
        station_name = 'PERTH';
    end
    [range_pred,Az_pred,El_pred] = ECI2RAzEl(station_name,et_vis_global(i),yy_sol_b(:,i));
    meas_pred = [range_pred; Az_pred; El_pred];
    meas_real = meas_real_global(:,i);
    residual_sol_b(:,i) = W_meas*(meas_pred - meas_real);    % residual = weight*(pred_meas - real_meas)
end
fprintf('Process took %f seconds\n',toc);
%% EX2.3c
% propagate refence state;
parameters_c = parameters_b;
parameters_c.motion = 'J2';

% propagate reference state
if et_vis_global(1) - t0 == 0
    [~,yy_ref] = ode113(@(t,y) eom_2body_J2(t,y,mu), et_vis_global - t0, state0, opt_ode);
    yy_ref_c = yy_ref(:,:)';
elseif et_vis_global(1) - t0 ~= 0
    [~,yy_ref] = ode113(@(t,y) eom_2body_J2(t,y,mu), [0, et_vis_global - t0], state0, opt_ode);
    yy_ref_c = yy_ref(2:end,:)';
end

% initialize residuals vector
fprintf('Point 2.3c: computing residuals with respect to reference trajectory.\n');
tic;
residual_ref_c = zeros(3,length(et_vis_global));   
% start a loop to compute residuals at each epoch:
for i = 1:length(et_vis_global)
    if station_ID_global(i) == 1
        station_name = 'KOUROU';
    elseif station_ID_global(i) == 2
        station_name = 'PERTH';
    end
    [range_pred,Az_pred,El_pred] = ECI2RAzEl(station_name,et_vis_global(i),yy_ref_c(:,i));
    meas_pred = [range_pred; Az_pred; El_pred];
    meas_real = meas_real_global(:,i);
    residual_ref_c(:,i) = W_meas*(meas_pred - meas_real);    % residual = weight*(pred_meas - real_meas)
end
fprintf('Process took %f seconds\n',toc);


% Solve nonlinear least squares problem
opt_lsq = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter-detailed');
fun_c = @(x) costfunction_radar(x,parameters_c);
%call solver
fprintf('Point 2.3c: optimization process using measurements at Perth and Kourou with J2 perturbation.\n'); 
tic;
[x_sol_c, resnorm_c, residual_c, flag_c, ~, ~, Jac_c] = lsqnonlin(fun_c,state0, [], [], opt_lsq);
fprintf('Optimization took %f seconds.\n',toc);
% estimation of covariance matrix:
Jac_c = full(Jac_c);
P_ls_c = resnorm_c/( length(residual_c) - length(state0) ) .*inv(Jac_c'*Jac_c);

%% EX2.3d J2 motion and Perth only measurements: not needed but used to justify conclusions
% initialize residuals vector
parameters_d = parameters_a;
parameters_d.motion = 'J2';

% propagate reference state
if Perth.et_vis(1) - t0 == 0
    [~,yy_ref] = ode113(@(t,y) eom_2body_J2(t,y,mu), Perth.et_vis - t0, state0, opt_ode);
    yy_ref_d = yy_ref(:,:)';
elseif Perth.et_vis(1) - t0 ~= 0
    [~,yy_ref] = ode113(@(t,y) eom_2body_J2(t,y,mu), [0, Perth.et_vis - t0], state0, opt_ode);
    yy_ref_d = yy_ref(2:end,:)';
end

% initialize residuals vector
fprintf('Point 2.3d: computing residuals with respect to reference trajectory.\n');
tic;
residual_ref_d = zeros(3,length(Perth.et_vis));   
% start a loop to compute residuals at each epoch:
for i = 1:length(Perth.et_vis)
    [range_pred,Az_pred,El_pred] = ECI2RAzEl('PERTH',Perth.et_vis(i),yy_ref_d(:,i));
    meas_pred = [range_pred; Az_pred; El_pred];
    meas_real = Perth.sgp4.measures_noise(:,i);
    residual_ref_d(:,i) = W_meas*(meas_pred - meas_real);    % residual = weight*(pred_meas - real_meas)
end
fprintf('Process took %f seconds\n',toc);


% Solve nonlinear least squares problem
opt_lsq = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt','Display','iter-detailed');
fun_c = @(x) costfunction_radar(x,parameters_d);
%call solver
fprintf('Point 2.3d: optimization process using measurements at Perth with J2 perturbation.\n'); 
tic;
[x_sol_d, resnorm_d, residual_d, flag_d, ~, ~, Jac_d] = lsqnonlin(fun_c,state0, [], [], opt_lsq);
fprintf('Optimization took %f seconds.\n',toc);
% estimation of covariance matrix:
Jac_d = full(Jac_d);
P_ls_d = resnorm_d/( length(residual_d) - length(state0) ) .*inv(Jac_d'*Jac_d);

%% PLOT RESULTS AS COVARIANCE ELLIPSOIDS
delta_sol_a = state0 - x_sol_a;
delta_sol_b = state0 - x_sol_b;
delta_sol_c = state0 - x_sol_c;
delta_sol_d = state0 - x_sol_d;
% resulting position vector and covariance ellipse in 2D
ellipse_ls_a_xy = Cov_ellipse2D(3,delta_sol_a(1:2),P_ls_a(1:2,1:2));
ellipse_ls_b_xy = Cov_ellipse2D(3,delta_sol_b(1:2),P_ls_b(1:2,1:2));
ellipse_ls_c_xy = Cov_ellipse2D(3,delta_sol_c(1:2),P_ls_c(1:2,1:2));
ellipse_ls_d_xy = Cov_ellipse2D(3,delta_sol_d(1:2),P_ls_d(1:2,1:2));
figure();
hold on;
grid on;
scatter(0,0,70,'o','filled','Color',"#0072BD");
h1 = scatter(delta_sol_a(1),delta_sol_a(2),40,[0.8500 0.3250 0.0980],'o','filled');
plot(ellipse_ls_a_xy(1,:),ellipse_ls_a_xy(2,:),'LineWidth',2,'Color',h1.CData);
h2 = scatter(delta_sol_b(1),delta_sol_b(2),40,[0.4660 0.6740 0.1880],'o','filled');
plot(ellipse_ls_b_xy(1,:),ellipse_ls_b_xy(2,:),'LineWidth',2,'Color',h2.CData);
h3 = scatter(delta_sol_c(1),delta_sol_c(2),40,'o','filled');
plot(ellipse_ls_c_xy(1,:),ellipse_ls_c_xy(2,:),'LineWidth',2,'Color',h3.CData);
h4 = scatter(delta_sol_d(1),delta_sol_d(2),40,[0.9290 0.6940 0.1250],'o','filled');
plot(ellipse_ls_d_xy(1,:),ellipse_ls_d_xy(2,:),'LineWidth',2,'Color',h4.CData);
legend('$\mathbf{x}_0^{ref}$','Case A','', ...
                              'Case B','', ...
                              'Case C','', ...
                              'Case D','', ...
                              'Location','best','FontSize',15);
xlabel('$\Delta x$ [km]'); ylabel('$\Delta y$ [km]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Least squares solutions','FontSize',15);

% resulting position vector and covariance ellipse in 2D
ellipse_ls_a_yz = Cov_ellipse2D(3,delta_sol_a(2:3),P_ls_a(2:3,2:3));
ellipse_ls_b_yz = Cov_ellipse2D(3,delta_sol_b(2:3),P_ls_b(2:3,2:3));
ellipse_ls_c_yz = Cov_ellipse2D(3,delta_sol_c(2:3),P_ls_c(2:3,2:3));
ellipse_ls_d_yz = Cov_ellipse2D(3,delta_sol_d(2:3),P_ls_d(2:3,2:3));
figure();
hold on;
grid on;
scatter(0,0,70,'o','filled','Color',"#0072BD");
h1 = scatter(delta_sol_a(2),delta_sol_a(3),40,[0.8500 0.3250 0.0980],'o','filled');
plot(ellipse_ls_a_yz(1,:),ellipse_ls_a_yz(2,:),'LineWidth',2,'Color',h1.CData);
h2 = scatter(delta_sol_b(2),delta_sol_b(3),40,[0.4660 0.6740 0.1880],'o','filled');
plot(ellipse_ls_b_yz(1,:),ellipse_ls_b_yz(2,:),'LineWidth',2,'Color',h2.CData);
h3 = scatter(delta_sol_c(2),delta_sol_c(3),40,'o','filled');
plot(ellipse_ls_c_yz(1,:),ellipse_ls_c_yz(2,:),'LineWidth',2,'Color',h3.CData);
h4 = scatter(delta_sol_d(2),delta_sol_d(3),40,[0.9290 0.6940 0.1250],'o','filled');
plot(ellipse_ls_d_yz(1,:),ellipse_ls_d_yz(2,:),'LineWidth',2,'Color',h4.CData);
legend('$\mathbf{x}_0^{ref}$','Case A','', ...
                              'Case B','', ...
                              'Case C','', ...
                              'Case D','', ...
    'Location','best','FontSize',15);
xlabel('$\Delta y$ [km]'); ylabel('$\Delta z$ [km]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Least squares solutions','FontSize',15);

% resulting position vector and covariance ellipse in 2D
ellipse_ls_a_xz = Cov_ellipse2D(3,delta_sol_a([1 3]),P_ls_a([1 3],[1 3]));
ellipse_ls_b_xz = Cov_ellipse2D(3,delta_sol_b([1 3]),P_ls_b([1 3],[1 3]));
ellipse_ls_c_xz = Cov_ellipse2D(3,delta_sol_c([1 3]),P_ls_c([1 3],[1 3]));
ellipse_ls_d_xz = Cov_ellipse2D(3,delta_sol_d([1 3]),P_ls_d([1 3],[1 3]));
figure();
hold on;
grid on;
scatter(0,0,70,'o','filled','Color',"#0072BD");
h1 = scatter(delta_sol_a(1),delta_sol_a(3),40,[0.8500 0.3250 0.0980],'o','filled');
plot(ellipse_ls_a_xz(1,:),ellipse_ls_a_xz(2,:),'LineWidth',2,'Color',h1.CData);
h2 = scatter(delta_sol_b(1),delta_sol_b(3),40,[0.4660 0.6740 0.1880],'o','filled');
plot(ellipse_ls_b_xz(1,:),ellipse_ls_b_xz(2,:),'LineWidth',2,'Color',h2.CData);
h3 = scatter(delta_sol_c(1),delta_sol_c(3),40,'o','filled');
plot(ellipse_ls_c_xz(1,:),ellipse_ls_c_xz(2,:),'LineWidth',2,'Color',h3.CData);
h4 = scatter(delta_sol_d(1),delta_sol_d(3),40,[0.9290 0.6940 0.1250],'o','filled');
plot(ellipse_ls_d_xz(1,:),ellipse_ls_d_xz(2,:),'LineWidth',2,'Color',h4.CData);
legend('$\mathbf{x}_0^{ref}$','Case a','', ...
                              'Case b','', ...
                              'Case c','', ...
                              'Case d','', ...
    'Location','best','FontSize',15);
xlabel('$\Delta x$ [km]'); ylabel('$ \Delta z$ [km]');
ax = gca; ax.XAxis.FontSize = 15; ax.YAxis.FontSize = 15;
title('Least squares solutions','FontSize',15);
%% PRINT RESULTS
Measurements = {'Perth';'All';'Perth';'All'};
Dynamics = {'Keplerian';'Keplerian';'J2 perturbed';'J2 perturbed'};
Residual = [resnorm_a; resnorm_b;resnorm_d;resnorm_c];
Delta    = [norm(state0 - x_sol_a); norm(state0 - x_sol_b); norm(state0 - x_sol_d); norm(state0 - x_sol_c)];
Tab = table(Measurements, Dynamics, Residual, Delta);
disp(Tab);


%% EX2.1 FUNCTIONS

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


function [rho,Az,El] = ECI2RAzEl(StationName,et_vec,x_vec)
% function that given an existing facility computes range, azimuth and elevation
% in its topocentric horizon reference frame from the s/c state in the 
% ECIJ2000 frame.
%
% INPUT:
% StationName = string containing the name of the observing station.
% et_vec      = [1xN] vector of epochs at which s/s state is computed.
% x_vec       = [6xN] matrix contaning the s/c states.
%
% OUTPUT:
% rho = [1xN] array containing the spacecraft radius coordinates [depends on input units].
% Az  = [1xN] array containing the spacecraft Azimuth coordinates [deg].
% El  = [1xN] array containing the spacecraft Elevation coordinates [deg].

% obtain azimuth and elevation for both stations at each propagated epoch:
Az = zeros(1,length(et_vec)); 
El = zeros(1,length(et_vec));
rho = zeros(1,length(et_vec)); 
TopoFrame = [StationName,'_TOPO'];
% start the loop:
for i = 1:length(et_vec)
    ECI2TOPO = cspice_sxform('J2000',TopoFrame, et_vec(i));
    % compute relative distance and velocity w.r.t. the station
    rv_station = cspice_spkezr(StationName, et_vec(i), 'J2000', 'NONE', 'EARTH');
    rv_station_sat_eci = x_vec(:,i) - rv_station;
    % convert into topocentric frame:
    rv_station_sat_topo = ECI2TOPO*rv_station_sat_eci;
    % compute azimuth and elevation
    [rho(i),Az(i),El(i)] = cspice_reclat(rv_station_sat_topo(1:3));
end

% now convert quantities in degrees:
Az = rad2deg(Az);
El = rad2deg(El);
end



function [et_vis,Az_vis,El_vis,i_visibility] = visibility_window(AzEl,et_vec,min_El)
% function that given a set of azimuth and elevation values determines when the spacecraft is
% visible, i.e its elevation is above a given treshold.
%
% INPUT: 
% AzEl   = [2xN] matrix of Azimuth and Elevation angles (first and second row respectively). [deg]
% et_vec = [1xN] array of epoch at which angles are provided.
% min_El = treshold under which the satellite is not visible from the station. [deg]
%
% OUTPUT
% et_vis   = [1xn] array containing epochs determining the beginning and ending of visibility.
% AzEl_vis = [1xn] matrix contaning angles during the visibility windows.
% windows_index = [mx1] array containing indexes associated to start and end of visibility windows.
%                       Indexes referred to whole epoch vector. Logical array.
% loop over Elevation values and set elements to zero if under the treshold
i_visibility = AzEl(2,:) > min_El;
et_vis = et_vec(i_visibility);
Az = AzEl(1,:);
El = AzEl(2,:);
El_vis = El(i_visibility);
Az_vis = Az(i_visibility);
% i_visibility = double(i_visibility);
% windows_index = 1:length(et_vec);
% windows_index = windows_index(ischange(i_visibility));
%windows_index = windows_index(windows_index~=0);
end


%% EX2.2 FUNCTIONS
function [r_eci, v_eci] = sgp4_propagator(sat_rec, et_vec, dPsi, dEps)
% function that propagates the s/c trajectory using sgp4 along the given epochs.
%
% INPUT:
% sat_rec = structure containing the s/c information obtained from a TLE.
% et_vec  = [1xN] array of epochs at which to compute the s/c state in ephemeris time.
% dPsi    = first reference value for the nutation correction.
% dEps    = second reference value for the precession and nutation correction.
%
% OUTPUT:
% r_eci  = [3xN] array containing the s/c position in ECIJ2000 frame at each epoch, [km].
% v_eci  = [3xN] array containing the s/c velocity in ECIJ2000 frame at each epoch, [km/s].

[year,mon,day,hr,min,sec] = invjday(sat_rec.jdsatepoch, sat_rec.jdsatepochf);
sgp4_epoch_str = sprintf('%d-%02d-%02d %02d:%02d:%02.6f', [year,mon,day,hr,min,sec]);
sgp4_epoch_et = cspice_str2et(sgp4_epoch_str);
sgp4_et_vec = (et_vec - sgp4_epoch_et)/60;

r_eci = zeros(3,length(et_vec));
v_eci = zeros(3,length(et_vec));

for i = 1:length(sgp4_et_vec)
    % call sgp4
    [~,rteme,vteme] = sgp4(sat_rec, sgp4_et_vec(i));
    % current epoch in centuries
    t_cent = cspice_unitim(et_vec(i), 'ET', 'TDT')/cspice_jyear()/100;
    % convert from teme to eci frame
    [r_eci(:,i), v_eci(:,i),~] = teme2eci(rteme, vteme, [0;0;0], t_cent, dPsi, dEps);
end

end


%% EX2.3 functions

function residual = costfunction_radar(state_ref, parameters)
% function that computes the residual between a nominal trajectory and a set
% of given measurements.
% We assume radar measurements, thus they will concern only the position.
%
% INPUT:
% state_ref                = [6x1] vector containing the reference time at initial epoch that we want to optimize.
% parameters               = struct having 8 fields:
% parameters.t_ref         = epoch associated to the reference state. [s]
% parameters.et_vis        = [1xN] vector containing the epochs at to propagate reference state. [s]
% parameters.station_ID    = [1xN_vis] array of station at which measure is computed: 1 for Kourou, 2
%                            for Perth.
% parameters.mu            = gravitational parameter of the assumed keplerian motion.
% parameters.measurements  = [3xN_vis] vector containing the available measurements. In our case it contains
%                           estimates of the spacecraft range, azimuth and elevation.
% parameters.motion        = string to define the motion of the spacecraft, either'Keplerian' or 'J2'.
% parameters.W_meas        = [3x3] weight matrix for measurements.
% OUTPUT:
% cost   = defined as the sum of the squares of each residual epsilon.
t_ref = parameters.t_ref;
et_vis = parameters.et_vis;
station_ID = parameters.station_ID;

opt_ode = odeset('RelTol',2.5e-14,'AbsTol',2.5e-14);
residual = zeros(3*length(et_vis),1);
switch parameters.motion
    case 'Keplerian'
        % propagate the actual reference state 
        if et_vis(1) - t_ref == 0
            [~,yy_ref] = ode113(@(t,y) eom_2body(t,y,parameters.mu), et_vis - t_ref, state_ref, opt_ode);
            yy_ref = yy_ref(:,:)';
        elseif et_vis(1) - t_ref ~= 0
            [~,yy_ref] = ode113(@(t,y) eom_2body(t,y,parameters.mu), [0, et_vis - t_ref], state_ref, opt_ode);
            yy_ref = yy_ref(2:end,:)';
        end
        % start a loop to compute residuals at each epoch:
        for i = 1:length(et_vis)
            % select station from which actual measurements are computed
            if station_ID(i) == 1
                station_name = 'KOUROU';
            elseif station_ID(i) == 2
                station_name = 'PERTH';
            end
            [range,Az,El] = ECI2RAzEl(station_name,et_vis(i),yy_ref(:,i));
            meas_pred = [range;Az;El];
            meas_real = parameters.measurements(:,i);
            residual(3*i-2:3*i) = parameters.W_meas*(meas_pred - meas_real);
        end
    case 'J2'
        % propagate the actual reference state 
        if et_vis(1) - t_ref == 0
            [~,yy_ref] = ode113(@(t,y) eom_2body_J2(t,y,parameters.mu), et_vis - t_ref, state_ref, opt_ode);
            yy_ref = yy_ref(:,:)';
        elseif et_vis(1) - t_ref ~= 0
            [~,yy_ref] = ode113(@(t,y) eom_2body_J2(t,y,parameters.mu), [0, et_vis - t_ref], state_ref, opt_ode);
            yy_ref = yy_ref(2:end,:)';
        end
        for i = 1:length(et_vis)
            % select station from which actual measurements are computed
            if station_ID(i) == 1
                station_name = 'KOUROU';
            elseif station_ID(i) == 2
                station_name = 'PERTH';
            end
            [range,Az,El] = ECI2RAzEl(station_name,et_vis(i),yy_ref(:,i));
            meas_pred = [range;Az;El];
            meas_real = parameters.measurements(:,i);
            residual(3*i-2:3*i) = parameters.W_meas*(meas_pred - meas_real);
        end
end

end

function [x_dot] = eom_2body_J2(t, x, mu)
% eom for restricted 2-body problem with the positional J2 perturbation
% INPUT:
%
% x   = state vector (r, v)                               [km],[km/s]
% t   = actual epoch in ephemeris time                    [s]
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
% position in ECI:
r_eci = x(1:3);
% conversion to ECEF:
rot = cspice_pxform('J2000','ITRF93',t);
r_ecef = rot*r_eci;
% compute J2 term:
J2 = 0.0010826269;
r_Earth = cspice_bodvrd('Earth','RADII',3);
r_Earth = r_Earth(1);
acc_J2_ecef = 3/2*mu*J2/r_norm^3*(r_Earth/r_norm)^2*( 5*(r_ecef(3)/r_norm)^2*r_ecef - [1;1;3].*r_ecef);
% transform it back to eci:
acc_J2_eci = rot'*acc_J2_ecef;
x_dot(4:6) = x_dot(4:6) + acc_J2_eci;
end


function [ellipse] = Cov_ellipse2D(N_std,mean,P)
% function that computes the 2D covariance ellipse with a given confidence interval
%
% INPUT:
% P = [2x2] covariance matrix
% N_std = number of standard deviations which defines the confidence interval
% mean  = mean value around which covariance ellipse is defined
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

