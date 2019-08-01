%% gen all
clear; close all; clc;

% gen_test_traj_PVTOL('unconu',100);
% gen_test_traj_PVTOL('uncon',100);
% gen_test_traj_PVTOL('ccm',100);

% gen_test_traj_PVTOL('unconu',250);
% gen_test_traj_PVTOL('uncon',250);
% gen_test_traj_PVTOL('ccm',250);

% gen_test_traj_PVTOL('unconu',500);
% gen_test_traj_PVTOL('uncon',500);
% gen_test_traj_PVTOL('ccm',500);
% 
% gen_test_traj_PVTOL('unconu',1000);
% gen_test_traj_PVTOL('uncon',1000);
gen_test_traj_PVTOL('ccm',1000);

%% sim all (TVLQR)

% test_dynamics_PVTOL_bulk('unconu','lqr',100);
% test_dynamics_PVTOL_bulk('uncon','lqr',100);
% test_dynamics_PVTOL_bulk('ccm','lqr',100);

% test_dynamics_PVTOL_bulk('unconu','lqr',250);
% test_dynamics_PVTOL_bulk('uncon','lqr',250);
% test_dynamics_PVTOL_bulk('ccm','lqr',250);

% test_dynamics_PVTOL_bulk('unconu','lqr',500);
% test_dynamics_PVTOL_bulk('uncon','lqr',500);
% test_dynamics_PVTOL_bulk('ccm','lqr',500);
% 
% test_dynamics_PVTOL_bulk('unconu','lqr',1000);
% test_dynamics_PVTOL_bulk('uncon','lqr',1000);
test_dynamics_PVTOL_bulk('ccm','lqr',1000);
    
%% sim all (MPC)

% test_dynamics_PVTOL_bulk('unconu','mpc',100);
% test_dynamics_PVTOL_bulk('uncon','mpc',100);
% test_dynamics_PVTOL_bulk('ccm','mpc',100);
% 
% test_dynamics_PVTOL_bulk('unconu','mpc',250);
% test_dynamics_PVTOL_bulk('uncon','mpc',250);
% test_dynamics_PVTOL_bulk('ccm','mpc',250);

% test_dynamics_PVTOL_bulk('unconu','mpc',500);
% test_dynamics_PVTOL_bulk('uncon','mpc',500);
% test_dynamics_PVTOL_bulk('ccm','mpc',500);
% 
% test_dynamics_PVTOL_bulk('unconu','mpc',1000);
% test_dynamics_PVTOL_bulk('uncon','mpc',1000);
% test_dynamics_PVTOL_bulk('ccm','mpc',1000);