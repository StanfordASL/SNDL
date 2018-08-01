%% gen all
clear; close all; clc;

% gen_test_traj_PVTOL(6,'unconu',250);
% gen_test_traj_PVTOL(6,'uncon',250);
% gen_test_traj_PVTOL(6,'ccm',250);
% 
% gen_test_traj_PVTOL(6,'unconu',500);
% gen_test_traj_PVTOL(6,'uncon',500);
% gen_test_traj_PVTOL(6,'ccm',500);
% 
% gen_test_traj_PVTOL(6,'unconu',1000);
% gen_test_traj_PVTOL(6,'uncon',1000);
% gen_test_traj_PVTOL(6,'ccm',1000);

%% sim all

% test_dynamics_PVTOL_bulk('unconu',0,250);
% test_dynamics_PVTOL_bulk('uncon',0,250);
% test_dynamics_PVTOL_bulk('ccm',0,250);
% 
% test_dynamics_PVTOL_bulk('unconu',0,500);
% test_dynamics_PVTOL_bulk('uncon',0,500);
% test_dynamics_PVTOL_bulk('ccm',0,500);
% 
% test_dynamics_PVTOL_bulk('unconu',0,1000);
% test_dynamics_PVTOL_bulk('uncon',0,1000);
% test_dynamics_PVTOL_bulk('ccm',0,1000);
%     
% test_dynamics_PVTOL_bulk('ccm',1,250);
% test_dynamics_PVTOL_bulk('ccm',1,500);
test_dynamics_PVTOL_bulk('ccm',1,1000);

