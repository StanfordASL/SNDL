clear; close all; clc;

%%%%% RUN sndl_startup first %%%%%

%% Input data

%X: N x n matrix of sample points (full constraint set, including demonstration tuples)
%U: N x m matrix of control inputs
%Y: N x n matrix of x_dot at sample points
%Xtr_i: N_tr x 1 vector of indices in X used for regression loss
%Xc_i: Nc x 1 vector of indices in X to enforce CCM constraints (dynamically updating)

load('PVTOL_data_train.mat'); %X,U,Y
load('PVTOL_data_val.mat'); %X_val, U_val, Y_val

%for regression
N_tr = 1000; %upto N_max (max number of demonstration tuples)
Xtr_i = round((linspace(1,N_max,N_tr))'); 

%for constraints (dynamically updating)
Nc = 200; %initial size of constraint set
Xc_i = randperm(size(X,1),Nc); Xc_i = Xc_i';

fprintf('N:%d, Nc = %d, N_tr = %d\n',size(X,1), length(Xc_i), N_tr);

%% Setup params

n = size(X,2);
m = size(U,2);
N = size(X,1);

load_PVTOL_params;

%% Call learner

learn_loop;
