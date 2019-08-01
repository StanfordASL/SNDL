function test_dynamics_PVTOL_bulk(model,ctrl,N_itr)

%Simulate model tracking for the generated trajectories

%%%%%%% Inputs %%%%%%%%
%model: 'uncon': unconstrained, regularized baseline model
%      'unconu': unconstrained, un-regularized baseline model
%      'ccm': CCM regularized model
%ctrl: 'lqr': TV-LQR controller
%      'mpc': MPC tracking
%      'ccm': CCM controller
%N_itr: #points used in training (demonstration) set for regression loss

%% Load trajectories

load(strcat('results_new/PVTOL_test_traj_',model,'_',num2str(N_itr),'_new','.mat'));
N_sim = length(Tref);

%% Load initial states

load('results/PVTOL_init_guess.mat','X0');

%% Initialize data struct

if strcmp(ctrl,'ccm') 
    DATA = cell(N_sim,5);
else
    DATA = cell(N_sim,4);
end

%% Simulate

do_plot = 0;

for idx = 1:N_sim
    disp(idx);
    
    if strcmp(ctrl,'lqr')
        
        [DATA{idx,1},DATA{idx,2},DATA{idx,3},DATA{idx,4}] = test_dynamics_PVTOL_LQR(N_itr,do_plot,model,Tref{idx},x_opt{idx},u_opt{idx},X0(idx,:)');
        
    elseif strcmp(ctrl,'mpc')
        
        [DATA{idx,1},DATA{idx,2},DATA{idx,3},DATA{idx,4}] = test_dynamics_PVTOL_MPC(N_itr,do_plot,model,Tref{idx},x_opt{idx},u_opt{idx},X0(idx,:)');
        
    else
        
        [DATA{idx,1},DATA{idx,2},DATA{idx,3},DATA{idx,4},DATA{idx,5}] = test_dynamics_PVTOL_CCM(N_itr,do_plot,Tref{idx},x_opt{idx},u_opt{idx},X0(idx,:)');
    end
    
end

%% Save

if strcmp(model,'uncon')
    PVTOL_UNCON_DATA = DATA;
    save(strcat('results/PVTOL_UNCON_',num2str(N_itr),'_',upper(ctrl),'.mat'),'PVTOL_UNCON_DATA');
    
elseif strcmp(model,'unconu')
    PVTOL_UNCONU_DATA = DATA;
    save(strcat('results/PVTOL_UNCONU_',num2str(N_itr),'_',upper(ctrl),'.mat'),'PVTOL_UNCONU_DATA');
    
elseif strcmp(model,'ccm')
    PVTOL_CCM_DATA = DATA;
    save(strcat('results_new/PVTOL_CCM_',num2str(N_itr),'_',upper(ctrl),'_new','.mat'),'PVTOL_CCM_DATA');
end

 
end

