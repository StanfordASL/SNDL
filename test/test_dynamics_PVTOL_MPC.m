function [t_ref,errors,t_span,X] = test_dynamics_PVTOL_MPC(N_itr,do_plot,model,t_ref,x_opt,u_opt,x0)

%Simulate with MPC controller

close all;

%% Load model

if strcmp(model,'ccm')
    load(strcat('PVTOL_Dyn_Functions_',num2str(N_itr),'.mat')); %f, B
else
    load(strcat('PVTOL_',model,'_Dyn_Functions_',num2str(N_itr),'.mat')); %f, B
end

%% constants

n = 6;
m = 2;

%% True dynamics

dyn_fnc = @pvtol_dyn;

%% Nominal trajectory: (x_opt, u_opt)

dT_0 = t_ref(2)-t_ref(1); %resolution of reference trajectory
dT = dT_0/2; %resolution of simulation
T = t_ref(end);
x_opt = interp1(t_ref,x_opt,0:dT:T);
u_opt = interp1(t_ref,u_opt,0:dT:T);
t_ref = 0:dT:T;

dt = dT/2; %resolution of result
t_span = 0:dt:T;

%% Setup Tracking MPC prob

x_con = [-inf, inf;
         -inf, inf;
         -70*pi/180, 70*(pi/180);
         -5, 5;
         -5, 5;
         -4*pi/6,4*pi/6];
        
u_con = [1, 10;
         1, 10];
     
T_mpc = 2;   %mpc lookahead   
N_mpc = round(T_mpc/dT); %number of points to use from nominal traj
N_traj = round(6*T_mpc); %number of collocation points for MPC prob
if (mod(N_traj,2)==1)
    N_traj = N_traj+1;
end

if (T > T_mpc) %only setup MPC if T > T_mpc
    
    [MPC_Prob,L_e_full,t_grid] = ...
        setup_MPC(n,m,...
        f_h,B_h,df_h,dB_h,x_con,u_con,...
        N_traj,T_mpc,dT,...
        2*eye(n),0.1^2,(1e-4)^2,...
        8.0*eye(n),x_opt(1:N_mpc+1,:),eye(m),u_opt(1:N_mpc+1,:),t_ref(1:N_mpc+1),'MPC');
    
    %define warm soln struct
    mpc_result = struct('sol',0);
    tic
    [~,~,mpc_solved,mpc_result] = ...
        compute_MPC(MPC_Prob,x0,n,m,N_traj,L_e_full,[x_opt(1:N_mpc+1,:),u_opt(1:N_mpc+1,:)],t_ref(1:N_mpc+1),t_grid,x_con,u_con,mpc_result);
    fprintf('MPC:%d, solve_time:%.4f\n',mpc_solved,toc);
    
else %use LQR
    mpc_solved = 0;
end

%if MPC failed on first time-step, compute LQR
if (mpc_solved~=1)
    L_lqr = solve_TVLQR(t_ref,x_opt,u_opt,f_h,B_h,df_h,dB_h,n,m);
end

delta = 0.5;
solve_MPC = (0:delta:T)';
if solve_MPC(end)==T
   T_steps_MPC = length(solve_MPC)-1;
else
   T_steps_MPC = length(solve_MPC);
end
        
%Store MPC solution
MPC_state = cell(T_steps_MPC,1);
MPC_ctrl = cell(T_steps_MPC,1);
X_MPC = nan(T_steps_MPC,n);

%% Simulate

X = zeros(length(t_span),n);
U = zeros(length(t_span)-1,m);
U_FB = zeros(length(t_ref)-1, m);

errors = zeros(length(t_ref)-1,n);

%initialize
X(1,:) = x0';
x = x0;

%ODE settings
ode_options = odeset('RelTol',1e-6,'AbsTol',1e-9);

%simulate
i_mpc = 0;
ii = 0;
for i = 1:length(t_ref)-1
    
    if (mod( round(t_ref(i)*10000)/10000 ,delta)==0)
        fprintf('t:%.3f/%.3f \n',t_ref(i),T);
    end
    
    if (mod( round(t_ref(i)*10000)/10000 ,delta)==0 && mpc_solved==1)
        %solve MPC only if previous solve did not fail
        %figure out correct nominal range
        if (i+N_mpc <= length(t_ref))
            i_end = i+N_mpc;
            XU_nom = [x_opt(i:i_end,:),u_opt(i:i_end,:)];
        else
            i_end = length(t_ref);
            XU_nom = [x_opt(i:i_end,:),u_opt(i:i_end,:);
                    repmat(x_opt(end,:),i+N_mpc-i_end,1),...
                    repmat(u_opt(end,:),i+N_mpc-i_end,1)];
        end
        
        tic
        [MP_state,MP_ctrl,mpc_solved,mpc_result] = ...
            compute_MPC(MPC_Prob,x,n,m,N_traj,L_e_full,XU_nom,t_ref(1:N_mpc+1),t_grid,x_con,u_con,mpc_result);
        fprintf('t: %.3f, MPC:%d, solve_time: %.4f\n',t_ref(i),mpc_solved,toc);
        i_mpc = i_mpc + 1;
        X_MPC(i_mpc,:) = x';
        
        if (mpc_solved==1)  
            
            %record solution
            MPC_state{i_mpc} = MP_state;
            MPC_ctrl{i_mpc} = MP_ctrl;
        
            %reset index for MPC control
            ii = 1;
        else
            %switch to LQR
            T_steps_MPC = i_mpc - 1;
            L_lqr = solve_TVLQR(t_ref(i:end)-t_ref(i),x_opt(i:end,:),u_opt(i:end,:),f_h,B_h,df_h,dB_h,n,m);
            ii = 1;
        end
    else
        ii = ii+1;
    end
    
    %nominal state and control
    x_nom = x_opt(i,:)';
    u_nom = u_opt(i:i+1,:); %[t_ref(i), t_ref(i+1)]
    
    %record error
    errors(i,:) = x_nom'-x';
    
    %compute feedback:
    if (mpc_solved==1)
        u_fb = MP_ctrl(ii,:)';
    else
        u_fb = L_lqr(:, :, ii)*(x - x_nom);
    end
    
    %store CCM control
    U_FB(i,:) = u_fb';
    
    %net control:
    U(1+(i-1)*(dT/dt):i*(dT/dt),:) = repmat(u_nom(1,:)+u_fb',(dT/dt),1);
    
    %simulate true dynamics
    [~, state_sim] = ode113(@(t_sim,state_sim) ode_sim(t_sim,state_sim,[t_ref(i),t_ref(i+1)],u_nom,u_fb,...
        dyn_fnc),[t_ref(i):dt:t_ref(i+1)],x,ode_options);
    
    %record
    state_sim(:,3) = wrapToPi(state_sim(:,3));
    X(1+(i-1)*(dT/dt):i*(dT/dt)+1,:) = state_sim;
    
    %update
    x = state_sim(end,:)';
    
end

%% Plots

if (do_plot)
    figure()
    hold on;
    plot(x_opt(:,1),x_opt(:,2),'r-','linewidth',2);
    plot(X(:,1),X(:,2),'b-','linewidth',2);
    plot(X_MPC(:,1),X_MPC(:,2),'ko','markersize',6,'markerfacecolor','k');
    plot(X(1,1),X(1,2),'go','markerfacecolor','g');
    plot(X(end,1),X(end,2),'ro','markerfacecolor','r');
    
    for i_mpc = 1:T_steps_MPC
        plot(MPC_state{i_mpc}(:,1),MPC_state{i_mpc}(:,2),'k-','linewidth',1.5);
        plot(MPC_state{i_mpc}(end,1),MPC_state{i_mpc}(end,2),'ko','markersize',6,'markerfacecolor','k');
    end
    
    for i = 1:(delta/dt):size(X,1)
        plot_pvtol(X(i,1),X(i,2),X(i,3),'b');
    end
    for i = 1:(delta/dT):size(x_opt,1)
        plot_pvtol(x_opt(i,1),x_opt(i,2),x_opt(i,3),'r');
    end
    grid on
    xlabel('X'); ylabel('Z');
    
    figure()
    subplot(2,1,1)
    plot(t_ref,u_opt,'linewidth',2);
    grid on
    xlabel('Time [s]'); title('Nominal control');
    
    subplot(2,1,2)
    plot(t_ref(1:end-1),U_FB,'linewidth',2);
    grid on
    xlabel('Time [s]'); title('Feedback control');
    
end

end

function dx_dt = ode_sim(t,x,dt_span,u_nom,u_fb,dyn_fnc)

u = interp1(dt_span,u_nom,t) + u_fb';

%clip control
for i = 1:2
    u(i) = max(min(u(i),10),1);
end

dx_dt = dyn_fnc(0,x,u');

end


