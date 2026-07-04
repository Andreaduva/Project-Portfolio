%% ex1 

clear;
close all;
clc;

set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% defining useful parameters
fun = @(x) cos(x) - x;
x_val = linspace(-2,2,100);
y_val = fun(x_val);
a = -5;
b = 5;
tol = 1e-8;
iter_max = 200;
% plot the function graph around the zero
figure();
p1 = plot(x_val,y_val,'b','LineWidth',5);
hold on;
xline(0,'LineWidth',2); yline(0,'LineWidth',2);
xlabel('$x$', 'FontSize',16);
ylabel('$y$', 'FontSize',16);
grid on;
legend(p1,'$f(x) = \cos(x)-x$','FontSize',14);
xlim([-2 2]);

for j = 1:7  % run it 7 times to get a stable computational time
tic;
[x_bis,iter_bis] = bisection_MSAS(fun,a,b,tol,iter_max);
t_bisec = toc;
tic;
[x_sec,iter_sec] = secant_MSAS(fun,a,b,tol,iter_max);
t_secant = toc;
tic;
[x_rf,iter_rf] = regula_falsi_MSAS(fun,a,b,tol,iter_max);
t_falsi = toc;
end

%% ex2

% contour plot to identify the zeros
f1 = @(x1,x2) x1.^2-x1-x2;
f2 = @(x1,x2) (x1.^2)/16 + x2.^2 -1;

x1 = linspace(-10,10,100);
x2 = linspace(-10,10,100);
[X,Y] = meshgrid(x1,x2);
figure();
contour(X,Y,f1(X,Y),[0 0],'r','ShowText','on','LineWidth',3);
hold on;
contour(X,Y,f2(X,Y),[0,0],'b','ShowText','on','LineWidth',3);
xline(0,'LineWidth',2); yline(0,'LineWidth',2);
xlabel('${x_1}$','FontSize',16);
ylabel('${x_2}$','FontSize',16);
legend('$f_1$','$f_2$','FontSize',14)
title('Contour lines for $f_1(x_1,x_2) = x_1^2-x_1-x_2 = 0$ and $f_2(x_1,x_2) = \frac{x_1^2}{16}+x_2^2-1=0$','FontSize',16);
grid on; xlim([-4 4]); ylim([-4 4]);


% initialize variables
x_guess1 = [3;3];
x_guess2 = [-0.5;1];
fun = @(x)  [x(1).^2-x(1)-x(2); x(1).^2/16 + x(2).^2-1];
dfun = @(x) [2*x(1) - 1, -1; x(1)/8, 2*x(2)];

% using the analytic Jacobian
[sol1,n_iter1] = Newton_exact(fun,dfun,x_guess1,1e-7,50);
[sol2,n_iter2] = Newton_exact(fun,dfun,x_guess2,1e-7,50);

% using forward differences
[sol_fd1,n_iter_fd1] = Newton_fd(fun,x_guess1,1e-7,50);
[sol_fd2,n_iter_fd2] = Newton_fd(fun,x_guess2,1e-7,50);

% using centered differences
[sol_cd1,n_iter_cd1] = Newton_cd(fun,x_guess1,1e-7,50);
[sol_cd2,n_iter_cd2] = Newton_cd(fun,x_guess2,1e-7,50);

%% ex3

% defining variables and functions
f=@(x,t) x-t.^2+1;
x_0 = 0.5;
sol_ex = @(t) t.^2 + 2*t + 1 - 0.5*exp(t);
t = [0 2];
h1 = 0.5;
h2 = 0.2;
h3 = 0.05;
h4 = 0.01;

% Heun's method
for i = 1:7               % 7 times to get a stable computational time  
    tic; [Heun_1,T_1] = Heun_MSAS(f,x_0,t,h1); t_Heun_1 = toc;
    tic; [Heun_2,T_2] = Heun_MSAS(f,x_0,t,h2); t_Heun_2 = toc;
    tic; [Heun_3,T_3] = Heun_MSAS(f,x_0,t,h3); t_Heun_3 = toc;
    tic; [Heun_4,T_4] = Heun_MSAS(f,x_0,t,h4); t_Heun_4 = toc;
end
% plot results
t_ex = linspace(0,2,50);
figure();
plot(t_ex,sol_ex(t_ex),'k','LineWidth',2);
hold on;
plot(T_1,Heun_1,'ob','MarkerSize',8,'LineWidth',2);
plot(T_2,Heun_2,'om','MarkerSize',8,'LineWidth',2);
plot(T_3,Heun_3,'or','MarkerSize',6,'LineWidth',1);
title('Heun''s\hspace{0.2cm}method','FontSize',18)
legend('Analytical', '$h_1 = 0.5$', '$h_2 = 0.2$', '$h_3 = 0.05$','FontSize',14,'Location','nw');
xlabel('$t$','interpreter','latex','FontSize',16);
ylabel('$x(t)$','interpreter','latex','FontSize',16);

%RK4 scheme
for i = 1:7
    tic; [RK4_1,~] = RK4_MSAS(f,x_0,t,h1); t_RK4_1 = toc;
    tic; [RK4_2,~] = RK4_MSAS(f,x_0,t,h2); t_RK4_2 = toc; 
    tic; [RK4_3,~] = RK4_MSAS(f,x_0,t,h3); t_RK4_3 = toc;
    tic; [RK4_4,~] = RK4_MSAS(f,x_0,t,h4); t_RK4_4 = toc;
end

% plot results
figure();
plot(t_ex,sol_ex(t_ex),'k','LineWidth',2);
hold on;
plot(T_1,RK4_1,'ob','MarkerSize',8,'LineWidth',2);
plot(T_2,RK4_2,'om','MarkerSize',8,'LineWidth',2);
plot(T_3,RK4_3,'or','MarkerSize',6,'LineWidth',1);
title('RK4\hspace{0.2cm}scheme','FontSize',16)
legend('Analytical', '$h_1 = 0.5$', '$h_2 = 0.2$', '$h_3 = 0.05$','FontSize',14,'Location','nw');
xlabel('$t$','FontSize',16);
ylabel('$x(t)$','FontSize',16);


% CPU time and integration error
int_er_H   = [ max(abs(Heun_1 - sol_ex(T_1))) , max(abs(Heun_2 - sol_ex(T_2))), max(abs(Heun_3 - sol_ex(T_3))), max(abs(Heun_4 - sol_ex(T_4))) ];
int_er_RK4 = [ max(abs(RK4_1 - sol_ex(T_1)))  , max(abs(RK4_2 - sol_ex(T_2))),  max(abs(RK4_3 - sol_ex(T_3))),  max(abs(RK4_4 - sol_ex(T_4)))  ];
CPU_time_H = [ t_Heun_1, t_Heun_2, t_Heun_3, t_Heun_4 ];
CPU_time_RK4 = [ t_RK4_1, t_RK4_2, t_RK4_3, t_RK4_4 ];
% plot results
figure();
hold on; 
plot(int_er_RK4,CPU_time_RK4,'ob',int_er_H,CPU_time_H,'or','MarkerSize',3,'LineWidth',3); grid on;
set(gca,'Xscale','log','Yscale','log');
xlabel('$max \mid x_{analytical}-x_{numeric} \mid$','FontSize',16); ylabel('CPU time$(s)$','FontSize',16);
title('CPU time $vs$ integration error','FontSize',20);
legend('RK4','Heun''s method', 'FontSize',16);

%% ex 4

% define the two operators:
A = @(alpha) [0 1; -1 2*cos(alpha)];
F_RK2 = @(h,alpha) eye(2) + h.*A(alpha) + 0.5*(h.*A(alpha))^2;
F_RK4 = @(h,alpha) eye(2) + h*A(alpha) + 0.5*(h*A(alpha))^2 + (1/6)*(h*A(alpha))^3 + (1/24)*(h*A(alpha))^4;
alpha_vect = linspace(pi,0,100);
% study the implicit function behaviour for alpha = [85°, 70°]
fun = @(h,alpha_i) max(abs(eig(F_RK4(h,alpha_i))))-1;
h = linspace(-5,5,50);
h_vect = zeros(1,50);
fun_vect_1 = zeros(1,50);
fun_vect_2 = zeros(1,50);
for j = 1:50
    fun_vect_1(j) = fun(h(j),deg2rad(85));
    fun_vect_2(j) = fun(h(j),deg2rad(70));
end
figure();
plot(h,fun_vect_1,h,fun_vect_2,'LineWidth',3);
xline(0,'LineWidth',2); yline(0,'LineWidth',2);
xlabel('$h$','FontSize',16); ylabel('$f(h)$','FontSize',16);
legend('$\alpha_1 = 85^{\circ}$','$\alpha_2 = 70^{\circ}$','FontSize',16);
title('Behaviour of the implicit function for two different $\alpha_{1,2}$','FontSize',20);
hold on;
ylim([-1 5]);
xlim([0 5]);

% call the solver
x0 = 2;
h_max_RK2 = h_alpha_solver(F_RK2, alpha_vect, x0);
x0 = 2.8;
h_max_RK4 = h_alpha_solver(F_RK4, alpha_vect, x0);

% we can derive one of the two complex conjugate eigenvetors of A
% analytically:
lambda_vect = cos(alpha_vect) + 1i*sin(alpha_vect);
% plot results
figure();
hold on;
plot(h_max_RK2.*lambda_vect,'r','LineWidth',3); plot(h_max_RK2.*conj(lambda_vect),'r','LineWidth',3);
plot(h_max_RK4.*lambda_vect,'b','LineWidth',3); plot(h_max_RK4.*conj(lambda_vect),'b','LineWidth',3);
%  better formulation of the problem
[X,Y] = meshgrid(linspace(-4,4,100));
z = X + 1i*Y;
R_4 = abs(1 + z + 0.5*(z).^2 + (1/6)*(z).^3 + (1/24)*(z).^4);
contour(X,Y,R_4,[1 1],'LineWidth',2);
% plot eigenvalues of problem 3
eig_3 = [0.5,0.2,0.05,0.01];
plot(eig_3,zeros(length(eig_3)),'xm','MarkerSize',16,'LineWidth',3);
grid on;
xline(0,'LineWidth',2); yline(0,'LineWidth',2);
xlim([-3.5 2]); ylim([-3 3]);
xlabel('Re$(h\lambda)$','FontSize',16); ylabel('Im$(h\lambda)$','FontSize',16);
title('Stability domain of Runge-Kutta schemes','FontSize',16');
legend('RK1','','RK2','','RK2 with different formulation','$h\lambda$ of exercise 3','FontSize',14);

%% ex 5

optRK1 = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt');
optRK2 = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt',...
    'functionTolerance',1e-12,'stepTolerance',1e-12);
optRK4 = optimoptions('lsqnonlin','Display','off','Algorithm','levenberg-marquardt',...
    'functionTolerance',1e-12,'stepTolerance',1e-12);

x0 = [1 ;1];
t = [0 1];
tol_vect = [1e-3, 1e-4, 1e-5, 1e-6];
alpha_vect = linspace(0,2*pi,350);
sol_an =@(A) expm(A)*x0;

h0RK1 = tol_vect;
h0RK2 = [0.05 0.015 4.7e-3,1.5e-3];
h0RK4 = [0.54,0.28,0.15,0.08];
h_RK1 = zeros(length(tol_vect),length(alpha_vect));
h_RK2 = zeros(length(tol_vect),length(alpha_vect));
h_RK4 = zeros(length(tol_vect),length(alpha_vect));

tic;
for i = 1:length(alpha_vect) % whole regions for alpha = [0, 2*pi]
    A = [ 0, 1; -1 2*cos(alpha_vect(i))];
    
    for j = 1:length(tol_vect)
    
    h_RK1(j,i) = lsqnonlin(@(h) norm(sol_an(A) - RK1_last(A, x0, t, h),Inf)-tol_vect(j),h0RK1(j),0,1,optRK1);
    h0RK1(j) = h_RK1(j,i);
    h_RK2(j,i) = lsqnonlin(@(h) norm(sol_an(A) - RK2_last(A, x0, t, h),Inf)-tol_vect(j),h0RK2(j),0,1,optRK2);
    h0RK2(j) = h_RK2(j,i);
    h_RK4(j,i) = lsqnonlin(@(h) norm(sol_an(A) - RK4_last(A, x0, t, h),Inf)-tol_vect(j),h0RK4(j),0,1,optRK4);
    h0RK4(j) = h_RK4(j,i);
    end
    
end
toc;

% extract values from alpha = pi
A = [ 0, 1; -1 2*cos(pi)];
h_RK1_pi = zeros(1,length(tol_vect));
h_RK2_pi = zeros(1,length(tol_vect));
h_RK4_pi = zeros(1,length(tol_vect));
h0RK1 = tol_vect;
h0RK2 = [0.05 0.015 4.7e-3,1.5e-3];
h0RK4 = [0.54,0.28,0.15,0.08];
for j = 1:length(tol_vect)   
h_RK1_pi(j) = lsqnonlin(@(h) norm(sol_an(A) - RK1_last(A, x0, t, h),Inf)-tol_vect(j),h0RK1(j),0,1,optRK1);
h_RK2_pi(j) = lsqnonlin(@(h) norm(sol_an(A) - RK2_last(A, x0, t, h),Inf)-tol_vect(j),h0RK2(j),0,1,optRK2);
h_RK4_pi(j) = lsqnonlin(@(h) norm(sol_an(A) - RK4_last(A, x0, t, h),Inf)-tol_vect(j),h0RK4(j),0,1,optRK4);
end


% plot result in h-lambda plane
lambda_vect = cos(alpha_vect) + 1i*sin(alpha_vect);
figure();
plot(h_RK1'.*lambda_vect');  % RK1
hold on; grid on;
xline(0); yline(0);
xlabel('Re$(h\lambda)$','FontSize',16); ylabel('Im$(h\lambda)$','FontSize',16);
legend('tol $= 10^{-3}$','tol $= 10^{-4}$','tol $= 10^{-5}$','tol $= 10^{-6}$','FontSize',16);
title('RK1 scheme','FontSize',18);

figure();
hold on;
plot(h_RK2'.*lambda_vect');  % RK2
xlabel('Re$(h\lambda)$','FontSize',16); ylabel('Im$(h\lambda)$','FontSize',16);
xline(0,'LineWidth',1); yline(0,'LineWidth',1);
grid on;
legend('tol $= 10^{-3}$','tol $= 10^{-4}$','tol $= 10^{-5}$','tol $= 10^{-6}$','FontSize',16);
title('RK2 scheme','FontSize',18);

figure();
hold on;
plot(h_RK4'.*lambda_vect');  % RK4
xlabel('Re$(h\lambda)$','FontSize',16); ylabel('Im$(h\lambda)$','FontSize',16);
xline(0,'LineWidth',1); yline(0,'LineWidth',1);
grid on;
legend('tol $= 10^{-3}$','tol $= 10^{-4}$','tol $= 10^{-5}$','tol $= 10^{-6}$','FontSize',16);
title('RK4 scheme','FontSize',18);

% function evaluation vs tol for alpha = pi
figure();
plot(ceil(1./h_RK1_pi),tol_vect,'ob', 2*ceil(1./h_RK2_pi),tol_vect,'dr',4*ceil(1./h_RK4_pi),tol_vect,'*m','MarkerSize',4,'LineWidth',3);
set(gca,'Xscale','log','Yscale','log'); grid on;
ylabel('tolerance','FontSize',16); xlabel('function evaluations','FontSize',16);
legend('RK1','RK2','RK4','FontSize',16);
title('function evaluations $vs$ tolerance for $\alpha = \pi$','FontSize',18);

%% ex6

% defining vectors and operators
alpha_vect = linspace(0,pi,100);
theta_vect = [0.1;0.3;0.4;0.7;0.9];
opt_solve = optimoptions('fsolve','Display','off');
A = @(alpha) [0 1; -1 2*cos(alpha)];
F_BI2_theta = @(h,alpha,theta) ( eye(2) -(1-theta)*h*A(alpha) + 0.5*(1-theta)^2*(h*A(alpha))^2 )\...
                               ( eye(2) +theta*h*A(alpha) + 0.5*(theta*h*A(alpha))^2 );
                           
h_sol = zeros(length(theta_vect),length(alpha_vect));
x0_vect = [2.5,5,10,5,2.5];

for j = 1:length(theta_vect)
    
    if theta_vect(j) < 0.6
        h_sol(j,:) = h_alpha_solver(@(h,alpha) F_BI2_theta(h,alpha,theta_vect(j)), alpha_vect, x0_vect(j));
    else
        h_sol(j,:) = h_alpha_solver(@(h,alpha) F_BI2_theta(h,alpha,theta_vect(j)), flip(alpha_vect), x0_vect(j));
        h_sol(j,:) = flip(h_sol(j,:));
    end
end

lambda_vect = cos(alpha_vect) + 1i*sin(alpha_vect);
figure();
hold on;
grid on;
p1 = plot(h_sol'.*lambda_vect');
colororder([p1(1).Color; p1(2).Color; p1(3).Color; p1(4).Color; p1(5).Color]);
p2 = plot(h_sol'.*conj(lambda_vect)');
xline(0,'LineWidth',1); yline(0,'LineWidth',1);
xlabel('Re$(h\lambda)$','FontSize',16); ylabel('Im$(h\lambda)$','FontSize',16);
title('Stability regions for $BI2_{\theta}$ schemes', 'FontSize',18);
legend('$\theta = 0.1$','$\theta = 0.3$','$\theta = 0.4$','$\theta = 0.7$','$\theta = 0.9$','FontSize',16);

%% ex 7

% define variables and functions
B = [ -180.5 219.5; 179.5 -220.5];
x0 = [1;1];
sol_an = @(t) expm(B*t)*x0;
state_fun = @(x,t) B*x;
h = 0.1;
theta = 0.1;
t = [0 5];
t_span = linspace(t(1),t(end),100);
x_an = zeros(2,length(t_span));
% analytical solution
for j = 1:length(t_span)
    x_an(:,j) = sol_an(t_span(j));
end
% numerical integration
[x_RK4,t_scheme] = RK4_MSAS(state_fun,x0,t,h);
[x_BI2,~]        = BI2_theta_lin(B,theta,x0,t,h);
% plot results
figure();
hold on;
grid on;
plot(t_span,x_an(1,:),'b','LineWidth',2);  % x1 behaviour
plot(t_scheme,x_BI2(1,:),'or','MarkerSize',6,'LineWidth',2);
xlabel('$t$','FontSize',16); ylabel('$x_1(t$)','FontSize',16);
legend('Analytical solution','$BI2_0.1$','FontSize',16);
title('Comparison between analytical and numerical solution for $x_1(t)$','FontSize',20);

figure();
hold on;
grid on;
plot(t_span,x_an(2,:),'b','LineWidth',2);   % x2 behaviour
plot(t_scheme,x_BI2(2,:),'or','MarkerSize',6,'LineWidth',2);
xlabel('$t$','FontSize',16); ylabel('$x_2(t$)','FontSize',16);
legend('Analytical solution','$BI2_0.1$','FontSize',16);
title('Comparison between analytical and numerical solution for $x_2(t)$','FontSize',20);

% stability domain and eigenvaues
figure();
hold on; grid on;
contour(X,Y,R_4,[1 1],'LineWidth',2);
p1 = plot(h_sol(1,:)'.*lambda_vect','m');
colororder([p1(1).Color]);
plot(h_sol(1,:)'.*conj(lambda_vect)','m');
plot( 0.1*eig(B),[ 0 ; 0],'xb','MarkerSize',6,'LineWidth',2);
xline(0,'LineWidth',1); yline(0,'LineWidth',1);
xlabel('Re$(h\lambda)$','FontSize',16); ylabel('Im$(h\lambda)$','FontSize',16);
title('Stability regions for $BI2_{0.1}$ and RK4 schemes', 'FontSize',18);
legend('RK4','$BI2_{0.1}$','','$h\lambda_{i}$','FontSize',16,'location','nw');
ylim([-3 3]); xlim([-40 3]);
%% functions ex1

function [zero,n_iter] = bisection_MSAS(fun,a,b,tol,n_max)

% security check
if fun(a)*fun(b) >= 0
    error('Provided a and b points need to satisfy: f(a)*f(b) < 0\n');
end

% define initial error and n_iteration to enter the method
err = tol + 1;
n_iter = 0;

%trial
f_a = fun(a);



while err > tol && n_iter < n_max
    x_mid = (a+b)/2;
    f_mid = fun(x_mid);
    if f_a*f_mid < 0
        b = x_mid;        
    else
        a = x_mid;
        f_a = f_mid;
    end
    n_iter = n_iter +1;
    err = 0.5*(b - a);
end

zero = x_mid; % solution is taken as the mid-point
      
end

function [zero, n_iter] = secant_MSAS(fun,x0,x1,tol,n_max)

% evaluation of the guesses:
f0 = fun(x0);
f1 = fun(x1);
% initializing error and number of iterations
err = tol + 1;
n_iter = 0;

while err > tol && n_iter < n_max
    
    x = x1 - (x1 - x0)/(f1 - f0)*f1;   % compute the new point
    
    err = abs(x1 - x);                 % update error and n_iteration
    n_iter = n_iter + 1;
    
    x0 = x1;                           % update values for next iteration
    x1 = x;
    f0 = f1;
    f1 = fun(x);
end

zero = x;

end

function [zero, n_iter] = regula_falsi_MSAS(fun,x0,x1,tol,n_max)


% for a better understanding it's better to rename the variables:
xm1 = x0;       % x at k-1.  At k+1 will be xp1
fm1 = fun(x0);  % corresponding function value
x0 = x1;        % x at k
f0 = fun(x1);   % corresponding function value

% security check
if fm1*f0 >= 0
    error('Provided x0 and x1 points need to satisfy: f(x0)*f(x1) < 0\n');
end

% define error and iteration number
err = tol + 1;
n_iter = 0;

while err > tol && n_iter < n_max
    
    xp1 = x0 - (x0 - xm1)/(f0 - fm1)*f0;   % x at k+1
    fp1 = fun(xp1);                        % function evaluation needed also for following comparison
    
    err = abs(fp1);                        % update error and n° of iterations
    n_iter = n_iter + 1;
    
    if fp1*f0 < 0                          % values just switch position by 1 and x at k-1 is released
        xm1 = x0;
        x0 = xp1;
        fm1 = f0;
        f0 = fp1;
    else                                   % the old value of x at k is released while x at k-1 is kept
        x0 = xp1;
        f0 = fp1;
    end


end
 zero = xp1;  % retrieve the solution
 
end


%% functions ex2 

function [zero,n_iter] = Newton_exact(fun,dfun,x0,tol,iter_max)

%define error, n° of iterations, and first x_k
err = tol+1;
n_iter = 0;
xk = x0;
fk = fun(xk);
    
while err > tol && n_iter < iter_max
        
    dfun_k = dfun(xk);     % evaluation of Jacobian matrix at x_k
    delta_xk = -dfun_k\fk; % evaluation of the step to take by solving a linear system instead of taking the inverse
    xkp1 = xk + delta_xk;  % compute x at k+1
        
    err = norm(delta_xk);  % update error and n° of iteration
    n_iter = n_iter + 1;
        
    xk = xkp1;             % update values for next iteration
    fk = fun(xk);
end
zero = xkp1;  % solution is the last iteration
    
end

function [zero_fd,n_iter_fd] = Newton_fd(fun,x0,tol,iter_max)

xk_fd = x0;                                 % initialize values
fk_fd = fun(xk_fd);

err_fd = tol + 1;
n_iter_fd = 0;                              % n° iterations for forward difference case

n_x = length(xk_fd);                           % how many indipendent variables enter fun()
dfun_k_fd = zeros(length(fk_fd),n_x);          % initialize Jacobian matrix
    
%forward differences    
while err_fd > tol && n_iter_fd < iter_max
        
    % approximating the Jacobian
    epsilon = max([sqrt(eps), sqrt(eps)*norm(xk_fd)]); % compute the perturbation       
    for i = 1:n_x
        x_pert_fd = xk_fd;                                 
        x_pert_fd(i) = xk_fd(i) + epsilon;                     % perturbed vector
        dfun_k_fd(:,i) = (1/epsilon)*(fun(x_pert_fd)-fun(xk_fd));   % forward difference
    end
    
    delta_xk_fd = -dfun_k_fd\fk_fd; % evaluation of the step to take
    xkp1_fd = xk_fd + delta_xk_fd;  % compute x at k+1
        
    err_fd = norm(delta_xk_fd);  % update error and n° of iteration
    n_iter_fd = n_iter_fd + 1;
        
    xk_fd = xkp1_fd;             % update values for next iteration
    fk_fd = fun(xk_fd);
                
end

zero_fd = xk_fd;              
 
end

function [zero_cd,n_iter_cd] = Newton_cd(fun,x0,tol,iter_max)
% centered differences
err_cd = tol + 1;
n_iter_cd = 0;

xk_cd = x0;
fk_cd = fun(xk_cd);

n_x = length(xk_cd);                           % how many indipendent variables enter fun()
dfun_k_cd = zeros(length(fk_cd),n_x);          % initialize Jacobian matrix

while err_cd > tol && n_iter_cd < iter_max
        
    % approximating the Jacobian
    epsilon = max([sqrt(eps), sqrt(eps)*norm(xk_cd)]); % compute the perturbation       
    for i = 1:n_x
        x_pertp_cd = xk_cd;                                 
        x_pertp_cd(i) = xk_cd(i) + epsilon;                     % perturbed vector x + eps
        x_pertm_cd = xk_cd;                                     % perturbed vector x - eps
        x_pertm_cd(i) = xk_cd(i) - epsilon; 
        dfun_k_cd(:,i) = (1/(2*epsilon))*(fun(x_pertp_cd)-fun(x_pertm_cd));   % centered difference
    end
    delta_xk_cd = -dfun_k_cd\fk_cd; % evaluation of the step to take
    xkp1_cd = xk_cd + delta_xk_cd;  % compute x at k+1
        
    err_cd = norm(delta_xk_cd);  % update error and n° of iteration
    n_iter_cd = n_iter_cd + 1;
        
    xk_cd = xkp1_cd;             % update values for next iteration
    fk_cd = fun(xk_cd);
                
end

zero_cd = xk_cd;                 % second column is centered differences result
end

%% numerical integrators

function [x_sol,t_sol] = Heun_MSAS(fun,x0,t,h)
% General-purpose fixed-step Heun's integrator.

% initialization of variables
t_sol = t(1):h:t(2);
n_step = length(t_sol);
x_sol = zeros(length(x0),n_step);
x_sol(:,1) = x0;
% enter the for loop
for i = 1:(n_step-1)
    tk = t_sol(i);
    xk = x_sol(:,i);
    fk1 = fun(xk,tk);                      % first function evaluation
    fk2 = fun(xk+h*fk1,tk+h);              % second function evaluation
    x_sol(:,i+1) = xk + 0.5*h*(fk1 + fk2);   % compute solution at next time step
end 
end

function [x_sol,t_sol] = RK4_MSAS(fun,x0,t,h)

% General-purpose fixed-step RK4 integrator.

% initialization of variables
t_sol = t(1):h:t(2);
n_step = length(t_sol);
x_sol = zeros(length(x0),n_step);
x_sol(:,1) = x0;
% enter the for loop
for i = 1:(n_step-1)
    tk = t_sol(:,i);
    xk = x_sol(:,i);
    fk1 = fun(xk,tk);                                         % first function evaluation
    fk2 = fun(xk + 0.5*h*fk1, tk + 0.5*h);                    % second function evaluation
    fk3 = fun(xk + 0.5*h*fk2, tk + 0.5*h);                    % third function evaluation
    fk4 = fun(xk +     h*fk3, tk +     h);                    % fourth function evaluation
    x_sol(:,i+1) = xk + (1/6)*h*(fk1 + 2*fk2 + 2*fk3 + fk4);      % compute solution at next time step
end
end

%% functions ex 5

function [x_RK1_last] = RK1_last(A,x0,t,h)
t_sol = t(1):h:t(2);
n_step = length(t_sol);
x_RK1_last = (eye(2)+h*A)^(n_step-1)*x0;
if t_sol(end)<t(end)
    h_end = t(end)-t_sol(end);
    x_RK1_last = (eye(2)+h_end*A)*x_RK1_last;
end

end

function [x_RK2_last] = RK2_last(A,x0,t,h)
t_sol = t(1):h:t(2);
n_step = length(t_sol);
x_RK2_last = (eye(2)+h*A+(h*A)^2/2)^(n_step-1)*x0;
if t_sol(end)<t(end)
    h_end = t(end)-t_sol(end);
    x_RK2_last = (eye(2)+h_end*A+(h_end*A)^2/2)*x_RK2_last;
end

end


function [x_RK4_last] = RK4_last(A,x0,t,h)
t_sol = t(1):h:t(2);
n_step = length(t_sol);
x_RK4_last = (eye(2)+h*A+(h*A)^2/2+(h*A)^3/6+(h*A)^4/24)^(n_step-1)*x0;
if t_sol(end)<t(end)
    h_end = t(end)-t_sol(end);
    x_RK4_last = (eye(2)+h_end*A+(h_end*A)^2/2+(h_end*A)^3/6+(h_end*A)^4/24)*x_RK4_last;
end
end

%% h,alpha solver (ex5,6,7)
function [h_max] = h_alpha_solver(F_oper, alpha,x0)
%
%
% Function that solves the problem "find h>=0 such that
% max(abs(eig(F)))=1".
%
% INPUT
% F_oper : operator F passed as anonymous function F_oper(h,alpha).
% alpha  : vector of angles for which we want to compute the solution.
% x0     : initial guess for the first value of alpha.
%
% OUTPUT
% h_max  : vector with the solutions for each given alpha
%


fun = @(h,alpha_i) max(abs(eig(F_oper(h,alpha_i))))-1;   % defining the implicit function
opt = optimoptions('fsolve','Display','off');

h_max = zeros(1,length(alpha));

for i = 1:length(alpha)
    h_max(i) = fsolve(@(h) fun(h,alpha(i)),x0,opt);       % solve the implicit function and store the value        
    x0 = h_max(i);

end
end

%% function ex7

function [x_sol,t_sol] = BI2_theta_lin(A,theta,x0,t,h)

% BI2 integration scheme using the theta provided by the user.
% Note that this scheme only works for 2x2 LTI systems,for other types of ODEs
% it should be modified.

% initialization of variables
t_sol = t(1):h:t(2);
n_step = length(t_sol);
x_sol = zeros(length(x0),n_step);
x_sol(:,1) = x0;
for i = 1:n_step-1
    xk = x_sol(:,i);
    x_sol(:,i+1) = ( eye(2) - (1-theta)*h*A + 0.5*(1-theta)^2*(h*A)^2 )\( eye(2) + theta*h*A + 0.5*(theta*h*A)^2 )*xk;
end
end