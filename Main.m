clear; close all; clc;

%%%%% RUN sndl_startup first %%%%%

%% Input data

%X: N x n matrix of sample points (full constraint set, including demonstration tuples)
%U: N_max x m matrix of control inputs
%Y: N_max x n matrix of x_dot at sample points

load('PVTOL_data_train.mat'); %X,U,Y
load('PVTOL_data_val.mat'); %X_val, U_val, Y_val

%for regression
N_tr = 1000; %upto N_max (max number of demonstration tuples)
Xtr_i = round((linspace(1,N_max,N_tr))'); 

%for constraints (dynamically updating)
Nc = 250;
Xc_i = randperm(size(X,1),Nc); Xc_i = Xc_i';

fprintf('N:%d, Nc = %d, N_tr = %d\n',size(X,1), length(Xc_i), N_tr);

%% Pre-load (as an example)

% load(sprintf('learned_functions/PVTOL_raw_%d',N_tr),'P_f','JP_f_all','om_f',...
%                                                     'P_b','om_b',...
%                                                     'P_Wp_all','P_Wnp_all','JP_Wp_all','om_w_p','om_w_np',...
%                                                     'Xtr_i');
% JP_f = JP_f_all(:,:,Xc_i,:);
% P_Wp = P_Wp_all(Xc_i,:);
% P_Wnp = P_Wnp_all(Xc_i,:);
% JP_Wp = JP_Wp_all(Xc_i,:,:);

%% Setup params

n = size(X,2);
m = size(U,2);
N = size(X,1);

X_dot = reshape(Y(Xtr_i,:)',N_tr*n,1);

load_PVTOL_params;

data = struct('X',X,'U',U,'X_dot',X_dot);

%% Call learner

learn_loop;