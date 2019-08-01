function [t_ref,errors,t_span,X] = test_dynamics_PVTOL_LQR(N_itr,do_plot,model,t_ref,x_opt,u_opt,x0)

%Simulate with TV-LQR controller

close all;

%% Load model

if strcmp(model,'ccm')
    load(strcat('PVTOL_Dyn_Functions_',num2str(N_itr),'_new','.mat')); %f, B
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

%% Setup controller

%solve TV-LQR
L_lqr = solve_TVLQR(t_ref,x_opt,u_opt,f_h,B_h,df_h,dB_h,n,m);

%% Simulate

X = zeros(length(t_span),n);
U = zeros(length(t_span)-1,m);
U_ILQR = zeros(length(t_ref)-1, m);

errors = zeros(length(t_ref)-1,n);

%initialize
X(1,:) = x0';
x = x0;

%ODE settings
ode_options = odeset('RelTol',1e-6,'AbsTol',1e-9);

%simulate

for i = 1:length(t_ref)-1
    
    %nominal state and control
    x_nom = x_opt(i,:)';
    u_nom = u_opt(i:i+1,:); %[t_ref(i), t_ref(i+1)]
    
    %record error
    errors(i,:) = x_nom'-x';
    
    %compute feedback:
    u_fb = L_lqr(:, :, i)*(x - x_nom);
    
    %store CCM control
    U_ILQR(i,:) = u_fb';
    
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
    plot(x_opt(:,1),x_opt(:,2),'r-','linewidth',2);
    hold on
    plot(X(:,1),X(:,2),'b-','linewidth',2);
    plot(X(1,1),X(1,2),'go','markerfacecolor','g');
    plot(X(end,1),X(end,2),'ro','markerfacecolor','r');
    for i = 1:(0.5/dt):size(X,1)
        plot_pvtol(X(i,1),X(i,2),X(i,3),'b');
    end
    for i = 1:(0.5/dT):size(x_opt,1)
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
    plot(t_ref(1:end-1),U_ILQR,'linewidth',2);
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