%
clear; close all; clc

%% Load learned model

load('PVTOL_H_Dyn_Functions_1000_new.mat');
% load('PVTOL_H_uncon_Dyn_Functions_1000_new.mat');      

%load('PVTOL_H_Dyn_Functions_150.mat');
%load('PVTOL_H_uncon_Dyn_Functions_150.mat');     


%% Setup traj optimizer

addpath('Traj_Opt');

x_con = [-inf, inf;
         -inf, inf;
         -70*pi/180, 70*(pi/180);
         -5, 5;
         -5, 5;
         -pi, pi];
        
u_con = [2, 20;
         -3*pi/4, 3*pi/4];
     
n = 6; m = 2;

%% Create reference trajectory

ry = 1.0;
rz = 0.7;

dt = (1/250);
T = 10;
om = 2*pi/T;

alpha = (om^2)/(4*pi);
T_break = om/(2*alpha);

x0_ramp = [0;-1.0;zeros(4,1)];

%% Solve ramp

time_ramp = (0:dt:T_break)';
theta_ramp = alpha*(time_ramp.^2);

py_ramp = ry - ry*cos(theta_ramp) + x0_ramp(1);
pz_ramp = -rz * sin(2*theta_ramp) + x0_ramp(2);

vy_ramp = ry*sin(theta_ramp)*2*alpha.*time_ramp;
vz_ramp = -rz*cos(2*theta_ramp)*2*(2*alpha).*time_ramp;

xT_ramp = [py_ramp(end);pz_ramp(end);0;vy_ramp(end);vz_ramp(end);0];

N_traj_break = round(10*T_break);
if (mod(N_traj_break,2)==1)
    N_traj_break = N_traj_break+1;
end

[MP_Prob,L_e_full,t_grid] = ...
    setup_MP_Ref(n,m,...
    f_h,B_h,df_h,dB_h,x_con,u_con,...
    N_traj_break,T_break,dt,...
    eye(n),0.0025,0.0025,...
    100*diag([1;1;zeros(4,1)]),[py_ramp,pz_ramp],xT_ramp,...
    0.01*eye(2),[10;0],'MP');

XU_guess = [py_ramp,pz_ramp,zeros(length(time_ramp),1),vy_ramp,vz_ramp,zeros(length(time_ramp),1),...
            10*ones(length(time_ramp),1),zeros(length(time_ramp),1)];
    
fprintf('Starting traj opt...');
[x_opt_ramp,u_opt_ramp,qual,Result_ramp] = compute_MP_ref(MP_Prob,x0_ramp,n,m,N_traj_break,L_e_full,XU_guess,...
                                                            time_ramp,t_grid,x_con);
fprintf('Done, converged: %d \n',qual);

%% Check soln

figure(1)
subplot(3,1,1)
plot(py_ramp,pz_ramp,'b--','linewidth',3); hold on
plot(x_opt_ramp(:,1),x_opt_ramp(:,2),'r-','linewidth',1);
set(gca,'Ydir','reverse')

% keyboard;

%% Create flat

x0_flat = x_opt_ramp(end,:)';

time_flat = (0:dt:T)';
theta_flat = om*time_flat+pi;

py_flat = ry - ry*cos(theta_flat) + x0_ramp(1);
pz_flat = -rz * sin(2*theta_flat) + x0_ramp(2);

%% Check

figure(2)
plot(py_ramp,pz_ramp,'b-','linewidth',3); hold on
plot(py_flat,pz_flat,'r-');
set(gca,'Ydir','reverse');

% keyboard;

%% Solve flat

N_traj_flat = round(10*T);
if (mod(N_traj_flat,2)==1)
    N_traj_flat = N_traj_flat+1;
end

[MP_Prob,L_e_full,t_grid] = ...
    setup_MP_Ref(n,m,...
    f_h,B_h,df_h,dB_h,x_con,u_con,...
    N_traj_flat,T,dt,...
    eye(n),0.0025,0.0025,...
    100*diag([1;1;zeros(4,1)]),[py_flat,pz_flat],x0_flat,...
    0.01*eye(2),[10;0],'MP');

XU_guess = [py_flat,pz_flat,zeros(length(time_flat),4),...
            10*ones(length(time_flat),1),zeros(length(time_flat),1)];
    
fprintf('Starting traj opt...');
[x_opt_flat,u_opt_flat,qual,Result_flat] = compute_MP_ref(MP_Prob,x0_flat,n,m,N_traj_flat,L_e_full,...
                                                           XU_guess,time_flat,t_grid,x_con);
fprintf('Done, converged: %d \n',qual);

%% Check 

figure(1)
subplot(3,1,2)
plot(py_flat,pz_flat,'b--','linewidth',3); hold on
plot(x_opt_flat(:,1),x_opt_flat(:,2),'r-','linewidth',1);
set(gca,'Ydir','reverse')

% keyboard;

%% Create slow-down

x0_slow = x_opt_flat(end,:)';

time_slow = time_ramp;
theta_slow = 2*pi-alpha*((time_slow-T_break).^2);

py_slow = ry - ry*cos(theta_slow) + x0_ramp(1);
pz_slow = -rz * sin(2*theta_slow) + x0_ramp(2);

%% Check

figure(3)
plot(py_ramp,pz_ramp,'b-','linewidth',3); hold on
plot(py_flat,pz_flat,'r-');
plot(py_slow,pz_slow,'g-','linewidth',3);
set(gca,'Ydir','reverse');

% keyboard;

%% Solve slow-down

N_traj_slow = round(10*T_break);
if (mod(N_traj_slow,2)==1)
    N_traj_slow = N_traj_slow+1;
end

[MP_Prob,L_e_full,t_grid] = ...
    setup_MP_Ref(n,m,...
    f_h,B_h,df_h,dB_h,x_con,u_con,...
    N_traj_slow,T_break,dt,...
    eye(n),0.0025,0.0025,...
    100*diag([1;1;zeros(4,1)]),[py_slow,pz_slow],x0_ramp,...
    0.01*eye(2),[10;0],'MP');

XU_guess = [py_slow,pz_slow,zeros(length(time_slow),4),...
            10*ones(length(time_slow),1),zeros(length(time_slow),1)];
    
fprintf('Starting traj opt...');
[x_opt_slow,u_opt_slow,qual,Result_slow] = compute_MP_ref(MP_Prob,x0_slow,n,m,N_traj_slow,L_e_full,...
                                                           XU_guess,time_slow,t_grid,x_con);
fprintf('Done, converged: %d \n',qual);

%% Plot

figure(1)
subplot(3,1,3)
plot(py_slow,pz_slow,'b--','linewidth',3); hold on
plot(x_opt_slow(:,1),x_opt_slow(:,2),'r-','linewidth',1);
set(gca,'Ydir','reverse')

% keyboard;

%% Plot state & control solutions

figure()
subplot(2,1,1)
plot(time_ramp,x_opt_ramp,'linewidth',2);
legend('y','z','\phi','vy','vz','phi dot');
subplot(2,1,2)
plot(time_ramp,u_opt_ramp,'linewidth',2);
legend('f','phi dot cmd');

figure()
subplot(2,1,1)
plot(time_flat,x_opt_flat,'linewidth',2);
legend('y','z','\phi','vy','vz','phi dot');
subplot(2,1,2)
plot(time_flat,u_opt_flat,'linewidth',2);
legend('f','phi dot cmd');

figure()
subplot(2,1,1)
plot(time_slow,x_opt_slow,'linewidth',2);
legend('y','z','\phi','vy','vz','phi dot');
subplot(2,1,2)
plot(time_slow,u_opt_slow,'linewidth',2);
legend('f','phi dot cmd');

%% LQR Gains

L_lqr_ramp = solve_TVLQR(time_ramp,x_opt_ramp,u_opt_ramp,f_h,B_h,df_h,dB_h,n,m);
L_lqr_ramp_vec = (reshape(L_lqr_ramp(:),[],length(time_ramp)))';

L_lqr_flat = solve_TVLQR(time_flat,x_opt_flat,u_opt_flat,f_h,B_h,df_h,dB_h,n,m);
L_lqr_flat_vec = (reshape(L_lqr_flat(:),[],length(time_flat)))';

L_lqr_slow = solve_TVLQR(time_slow,x_opt_slow,u_opt_slow,f_h,B_h,df_h,dB_h,n,m);
L_lqr_slow_vec = (reshape(L_lqr_slow(:),[],length(time_slow)))';

%verify solution

figure()
subplot(3,1,1)
plot(time_ramp, L_lqr_ramp_vec);

subplot(3,1,2)
plot(time_flat, L_lqr_flat_vec);

subplot(3,1,3)
plot(time_slow, L_lqr_slow_vec);

save('FIG8_CCM_1000.mat');

keyboard;


%% Printing

C_r = zeros(N_traj_break+1,n);
for i = 1:n
    c_r = Result_ramp.x_k(i:n:n*(N_traj_break)+i);
    C_r(:,i) = c_r;
end

for i = 1:size(C_r,1)
    fprintf('%.8f,',C_r(i,:)); fprintf('\n');
end

fprintf('******************************\n');

C_f = zeros(N_traj_break+1,n);
for i = 1:n
    c_f = Result_flat.x_k(i:n:n*(N_traj_flat)+i);
    C_f(:,i) = c_f;
end

for i = 1:size(C_f,1)
    fprintf('%.8f,',C_f(i,:)); fprintf('\n');
end

fprintf('******************************\n');

C_s = zeros(N_traj_break+1,n);
for i = 1:n
    c_s = Result_slow.x_k(i:n:n*(N_traj_break)+i);
    C_s(:,i) = c_s;
end

for i = 1:size(C_s,1)
    fprintf('%.8f,',C_s(i,:)); fprintf('\n');
end

fprintf('******************************\n');

keyboard; 

%% print U

clc;

U_r = zeros(N_traj_break+1,m);
for j = 1:m
    u_r = Result_ramp.x_k(n*(N_traj_break+1)+j:m:end-(m-j));
    U_r(:,j) = u_r;
end

for i = 1:size(U_r,1)
    fprintf('%.8f,',U_r(i,:)); fprintf('\n');
end

fprintf('******************************\n');

U_f = zeros(N_traj_flat+1,m);
for j = 1:m
    u_f = Result_flat.x_k(n*(N_traj_break+1)+j:m:end-(m-j));
    U_f(:,j) = u_f;
end

for i = 1:size(U_f,1)
    fprintf('%.8f,',U_f(i,:)); fprintf('\n');
end

fprintf('******************************\n');

U_s = zeros(N_traj_break+1,m);
for j = 1:m
    u_s = Result_slow.x_k(n*(N_traj_break+1)+j:m:end-(m-j));
    U_s(:,j) = u_s;
end

for i = 1:size(U_s,1)
    fprintf('%.8f,',U_s(i,:)); fprintf('\n');
end



%% Convert Cheby approximations

cheby_order = 3;
cheby_pnts = 500;

C_lqr_ramp = zeros(6,cheby_order+1,2);
C_lqr_flat = zeros(6,cheby_order+1,2);
C_lqr_slow = zeros(6,cheby_order+1,2);

fprintf('getting cheby coeffs...');

for j = 1:2
    for i = 1:6
        fnc = @(t) interp1(time_ramp,reshape(L_lqr_ramp(j,i,:),[],1),t);
        C_lqr_ramp(i,:,j) = get_cheby_approx(fnc,[0,T_break],cheby_pnts,cheby_order);
        
        fnc = @(t) interp1(time_flat,reshape(L_lqr_flat(j,i,:),[],1),t);
        C_lqr_flat(i,:,j) = get_cheby_approx(fnc,[0,T],cheby_pnts,cheby_order);
        
        fnc = @(t) interp1(time_slow,reshape(L_lqr_slow(j,i,:),[],1),t);
        C_lqr_slow(i,:,j) = get_cheby_approx(fnc,[0,T_break],cheby_pnts,cheby_order);
    end
end
fprintf('done\n');

%% print TV-LQR cheby coeffs

clc;

for i = 1:6
    fprintf('%.6f,',C_lqr_ramp(i,:,1)); fprintf('\n');
end
fprintf('--------\n');
for i = 1:6
    fprintf('%.6f,',C_lqr_ramp(i,:,2)); fprintf('\n');
end
fprintf('****************\n');

for i = 1:6
    fprintf('%.6f,',C_lqr_flat(i,:,1)); fprintf('\n');
end
fprintf('--------\n');
for i = 1:6
    fprintf('%.6f,',C_lqr_flat(i,:,2)); fprintf('\n');
end
fprintf('****************\n');

for i = 1:6
    fprintf('%.6f,',C_lqr_slow(i,:,1)); fprintf('\n');
end
fprintf('--------\n');
for i = 1:6
    fprintf('%.6f,',C_lqr_slow(i,:,2)); fprintf('\n');
end
fprintf('****************\n');

%% Print static-LQR coeffs

clc

for j = 1:2
    fprintf('%.6f,', L_lqr_ramp(j,:,1)); fprintf('\n');
end
fprintf('********************\n');

for j = 1:2
    fprintf('%.6f,', L_lqr_flat(j,:,1)); fprintf('\n');
end
fprintf('********************\n');

for j = 1:2
    fprintf('%.6f,', L_lqr_slow(j,:,1)); fprintf('\n');
end
fprintf('********************\n');
