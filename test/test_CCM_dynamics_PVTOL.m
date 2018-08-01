function [t_ref,errors,t_span,X,E_geod,lamb] = test_CCM_dynamics_PVTOL(N_itr,do_plot,use_ccm,t_ref,x_opt,u_opt)

% Simulate CCM Model
close all; 

%% Load dynamics

load(strcat('PVTOL_Dyn_Functions_',num2str(N_itr),'.mat')); %f, B
load(strcat('PVTOL_Metric_Functions_',num2str(N_itr),'.mat')); %W, dW, lambda

%% constants

n = 6;
m = 2;
lambda = 0.1; %self chosen value - validate later.

%% True dynamics

dyn_fnc = @pvtol_dyn;
% dyn_fnc = @(t,x,u) f_h(x) + B_h(x)*u; %for sanity check only
    
%% Nominal trajectory: (x_opt, u_opt)

dT_0 = t_ref(2)-t_ref(1); %resolution of reference trajectory
dT = dT_0/2; %resolution of simulation
T = t_ref(end);
x_opt = interp1(t_ref,x_opt,0:dT:T);
u_opt = interp1(t_ref,u_opt,0:dT:T);
t_ref = 0:dT:T;

dt = dT/2; %resolution of result
t_span = 0:dt:T;

%% Setup CCM controller

if (use_ccm)
    %Setup CCM controller
    E_thresh = 0.001;
    x_lim_i = 3;
    x_lim = 75*(pi/180);
    geodesic_N = 3;
    [geo_Prob,cntrl_info] = setup_CCM_controller(n,m,geodesic_N,...
        x_lim_i,x_lim,...
        W_h,dW_h,f_h,B_h,df_h,...
        lambda,x_opt(1,:)',u_opt(1,:)',dT,E_thresh);
    ok_prev = 0; lambda_prev = lambda;
end

%% Setup sim

X = zeros(length(t_span),n);
U = zeros(length(t_span)-1,m);
U_FB = zeros(length(t_ref)-1,m);

E_geod = zeros(length(t_ref)-1,1);
lamb = zeros(length(t_ref)-1,1);
CCM_valid = zeros(length(t_ref)-1,1);
errors = zeros(length(t_ref)-1,n);

CCM_solved = -1*ones(length(t_ref)-1,1);
CCM_comp_time = zeros(length(t_ref)-1,1);

%initialize
X(1,:) = x_opt(1,:);
x = x_opt(1,:)';

%ODE settings
ode_options = odeset('RelTol',1e-6,'AbsTol',1e-9);

%pre-compute backup LQR controller
L_lqr = solve_TVLQR(t_ref,x_opt,u_opt,f_h,B_h,df_h,dB_h,n,m);

d_prev = 0;

%simulate
for i = 1:length(t_ref)-1
    
   %nominal state and control
   x_nom = x_opt(i,:)';
   u_nom = u_opt(i:i+1,:); %[t_ref(i), t_ref(i+1)]
   
   %record error
   errors(i,:) = x_nom'-x';
   
   %compute feedback:
   if (use_ccm)
       tic
       if sum(max(abs(x(x_lim_i)),abs(x_nom(x_lim_i)))<=x_lim)==length(x_lim_i)
           [E_geod(i),CCM_solved(i),geo_Prob,cntrl_info,u_fb,CCM_valid(i),lamb(i),d]  = compute_CCM_controller(geo_Prob,cntrl_info,x_nom,u_nom(1,:)',x,ok_prev,lambda_prev,d_prev);
       end
       CCM_comp_time(i) = toc;
       ok_prev = CCM_valid(i); lambda_prev = lamb(i); d_prev = d;
   end
   
   %if ccm not valid, switch to LQR
   if (~CCM_valid(i) || E_geod(i)<E_thresh)
       u_fb = L_lqr(:, :, i)*(x - x_nom);
   end
   
   %store FB control
   U_FB(i,:) = u_fb';
   
   %net control:
   U(1+(i-1)*(dT/dt):i*(dT/dt),:) = repmat( max(min(u_nom(1,:)+u_fb',10),1) ,(dT/dt),1);
   
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
    for i = 1:(0.5/dt):size(X,1)
        plot_pvtol(X(i,1),X(i,2),X(i,3),'b');
    end
    for i = 1:(0.5/dT):size(x_opt,1)
        plot_pvtol(x_opt(i,1),x_opt(i,2),x_opt(i,3),'r'); 
    end
    grid on
    plot(x_opt(1,1),x_opt(1,2),'go','markersize',10,'markerfacecolor','g');
    plot(x_opt(end,1),x_opt(end,2),'ro','markersize',10,'markerfacecolor','r');
    xlabel('X'); ylabel('Z');
    axis tight; axis equal;
    
    figure()
    subplot(2,1,1)
    plot(t_ref,u_opt,'linewidth',2);
    grid on
    xlabel('Time [s]'); title('Nominal control');
    
    subplot(2,1,2)
    plot(t_span(1:end-1),U,'linewidth',2);
    grid on
    xlabel('Time [s]'); title('Net control');
    
    figure()
    subplot(2,1,1)
    plot(t_ref(1:end-1),CCM_solved,'ro','markerfacecolor','b','markersize',10);
    grid on
    xlabel('Time [s]');
    title('solved');
    
    subplot(2,1,2)
    plot(t_ref(1:end-1),CCM_comp_time,'ro','markerfacecolor','b','markersize',10);
    grid on
    xlabel('Time [s]');
    title('solve time');
    
    figure()
    plot(t_ref(1:end-1),E_geod,'b-','linewidth',2);
    grid on
    xlabel('Time [s]');
    ylabel('Geodesic energy');
    
    figure()
    plot(t_ref(1:end-1),lamb,'b-','linewidth',2);
    hold on
    plot(t_ref(1:end-1),CCM_valid,'r-','linewidth',2);
    grid on
    xlabel('Time [s]');
    h_l = legend('$\lambda$','CCM valid');
    set(h_l,'interpreter','latex');
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