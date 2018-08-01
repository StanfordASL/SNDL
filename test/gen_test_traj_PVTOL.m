function gen_test_traj_PVTOL(N,type,N_itr)

%Generate test trajectories

%%%%%% Inputs %%%%%%%
%N: grid discretization over (px,pz) -> generate N^2 trajectories
%type: 'uncon': unconstrained, regularized baseline model
%      'unconu': unconstrained, un-regularized baseline model
%      'ccm': CCM regularized model
%N_itr: #points used in training (demonstration) set for regression loss

%% Load model

if strcmp(type,'uncon')
    load (strcat('PVTOL_uncon_Dyn_Functions_',num2str(N_itr),'.mat'));
elseif strcmp(type,'unconu')
    load (strcat('PVTOL_unconu_Dyn_Functions_',num2str(N_itr),'.mat'));
elseif strcmp(type,'ccm')
    load(strcat('PVTOL_Dyn_Functions_',num2str(N_itr),'.mat')); %f, B
else
    disp('unrecognized file type');
end

%% Setup traj optimizer

addpath('Traj_Opt');

state_con = [-inf, inf;
             -inf, inf;
            -70*pi/180, 70*(pi/180);
            -5, 5;
            -5, 5;
            -4*pi/6,4*pi/6];
        
u_con = [1, 10;
         1, 10];
     
x_eq = [0;0;0;0;0;0];  
u_eq = [0;0];

n = 6; m = 2;

%% Check for existing initial guesses

load_init = 0;
if exist('results/PVTOL_init_guess.mat','file')==2
    load('results/PVTOL_init_guess.mat','guess','X0');
    load_init = 1;
end

%% Gen

if (~load_init)
    x_range = linspace(-8,8,N);
    y_range = x_range;
    x0_range_np = [5*pi/180;3.5;3.5;10*pi/180]; %other than position states
    X0 = zeros(N^2,6);
    guess = cell(N^2,2);
end

Tref = cell(N^2,1);
x_opt = cell(N^2,1);
u_opt = cell(N^2,1);
qual = zeros(N^2,1);

DT = 0.01;

idx = 1;
for i = 1:N
    for j = 1:N
        
        if (~load_init)
            %generate initial start
            x0_p = [x_range(j);y_range(i)];
            x0_np = -x0_range_np + 2*x0_range_np.*rand(4,1);
            x0 = [x0_p;x0_np];
            X0(idx,:) = x0';
            
            %create initial spline guess
            [pos_coef,vel_coef,acc_coef,jer_coef,T_seg,coeff] = smooth_2D_path([x0_p';0,0],1,x0_np(2:3)');
            
            %simulate
            [Xs, Us, ~] =  simulate_PVTOL(pos_coef,vel_coef,acc_coef,jer_coef,T_seg{1},coeff{1});
            T = size(Xs,1)*0.1; XU_guess = [Xs,Us];
            XU_guess(end+1,:)=XU_guess(end,:);
            
            %store guess
            t_ref = (0:DT:T)';
            guess{idx,1} = t_ref; 
            guess{idx,2} = XU_guess;
            
        else
            %recall guess
            x0 = X0(idx,:)';
            t_ref = guess{idx,1};
            XU_guess = guess{idx,2};
        end
        
        %Setup and call traj optimizer
        N_traj = round(6*t_ref(end));
        if (mod(N_traj,2)==1)
            N_traj = N_traj+1;
        end
        [MP_Prob,L_e_full,t_grid] = ...
            setup_MP(n,m,...
               f_h,B_h(randn(n,1)),df_h,state_con,u_con,...
               N_traj,t_ref(end),DT,...
               2*eye(n),2e-1,1e-1,...
               zeros(n),x_eq,eye(2),u_eq,'MP');

        fprintf('Starting traj opt...');
        Tref{idx} = t_ref;
        [x_opt{idx},u_opt{idx},qual(idx)] = ...
            compute_MP(MP_Prob,x0,n,m,N_traj,L_e_full,XU_guess,0:0.1:t_ref(end),t_grid,state_con);
        fprintf('Done, converged: %d \n',qual(idx));
        %update count
        idx = idx+1;
    end
end

%save initial guesses
save('results/PVTOL_init_guess.mat','guess','X0');

%save all solutions
save(strcat('results/PVTOL_test_traj_',type,'_',num2str(N_itr),'.mat'),'Tref','x_opt','u_opt','qual');

%% Plot all traj
close all; 

figure()
hold all
for i = 1:N^2
    plot(x_opt{i}(:,1),x_opt{i}(:,2),'linewidth',2);
end
grid on

figure()
subplot(2,2,1)
hold all
for i = 1:N^2
    plot(x_opt{i}(:,3),'linewidth',2);
end
grid on
title('\theta');

subplot(2,2,2)
hold all
for i = 1:N^2
    plot(x_opt{i}(:,4),'linewidth',2);
end
grid on
title('v_x');

subplot(2,2,3)
hold all
for i = 1:N^2
    plot(x_opt{i}(:,5),'linewidth',2);
end
grid on
title('v_z');

subplot(2,2,4)
hold all
for i = 1:N^2
    plot(x_opt{i}(:,6),'linewidth',2);
end
grid on
title('\omega');

figure()
subplot(2,1,1)
plot(qual(:,1),'rx');
grid on





