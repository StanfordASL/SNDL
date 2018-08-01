function test_dynamics_PVTOL_bulk(type,hybrid,N_itr)

%Simulate model tracking for the generated trajectories

%%%%%%% Inputs %%%%%%%%
%type: 'uncon': unconstrained, regularized baseline model
%      'unconu': unconstrained, un-regularized baseline model
%      'ccm': CCM regularized model
%hybrid: 1 if want CCM-LQR hybrid controller (used only for CCM model sims)
%        0 if want LQR controller (used only for CCM model simulations)
%N_itr: #points used in training (demonstration) set for regression loss

%% Load trajectories

load(strcat('results/PVTOL_test_traj_',type,'_',num2str(N_itr),'.mat'));

%% Sim CCM
if strcmp(type,'ccm')
    %% ccm
    
    % setup data struct
    N_sim = length(Tref);
    PVTOL_CCM_DATA = cell(N_sim,6);
    
    % Go through and sim
    
    for idx = 1:N_sim
        disp(idx);
        [PVTOL_CCM_DATA{idx,1},...
            PVTOL_CCM_DATA{idx,2},...
            PVTOL_CCM_DATA{idx,3},...
            PVTOL_CCM_DATA{idx,4},...
            PVTOL_CCM_DATA{idx,5},...
            PVTOL_CCM_DATA{idx,6}] = test_CCM_dynamics_PVTOL(N_itr,0,hybrid,Tref{idx},x_opt{idx},u_opt{idx});
    end
    
    % save
    save(strcat('results/PVTOL_CCM_SIM_',num2str(N_itr),'_',num2str(hybrid),'.mat'),'PVTOL_CCM_DATA');
    
elseif strcmp(type,'uncon')
    %% uncon
    
    %setup data
    N_sim = length(Tref);
    PVTOL_UNCON_DATA = cell(N_sim,4);
    
    % Go through and sim
    
    for idx = 1:N_sim
        disp(idx);
        [PVTOL_UNCON_DATA{idx,1},...
            PVTOL_UNCON_DATA{idx,2},...
            PVTOL_UNCON_DATA{idx,3},...
            PVTOL_UNCON_DATA{idx,4}] = test_uncon_dynamics_PVTOL(N_itr,0,Tref{idx},x_opt{idx},u_opt{idx});
    end
    
    % save
    
    save(strcat('results/PVTOL_UNCON_SIM_',num2str(N_itr),'.mat'),'PVTOL_UNCON_DATA');
    
elseif strcmp(type,'unconu')
    %% uncon
    
    %setup data
    N_sim = length(Tref);
    PVTOL_UNCONU_DATA = cell(N_sim,4);
    
    % Go through and sim
    
    for idx = 1:N_sim
        disp(idx);
        [PVTOL_UNCONU_DATA{idx,1},...
            PVTOL_UNCONU_DATA{idx,2},...
            PVTOL_UNCONU_DATA{idx,3},...
            PVTOL_UNCONU_DATA{idx,4}] = test_unconu_dynamics_PVTOL(N_itr,0,Tref{idx},x_opt{idx},u_opt{idx});
    end
    
    % save
    
    save(strcat('results/PVTOL_UNCONU_SIM_',num2str(N_itr),'.mat'),'PVTOL_UNCONU_DATA');
    
end

