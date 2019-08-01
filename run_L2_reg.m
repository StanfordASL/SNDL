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

%% Setup params

n = size(X,2);
m = size(U,2);
N = size(X,1);

N_tr_list = [100,250,500,1000];

%% Do L2 check

val_reg_err_all = [];

for nn = 1:4
    N_tr = N_tr_list(nn);
    
    %% Load Features

    load(sprintf('learned_functions/PVTOL_raw_%d',N_tr),'P_f','JP_f_all','om_f',...
                                                        'P_b','om_b',...
                                                        'P_Wp_all','P_Wnp_all','JP_Wp_all','om_w_p','om_w_np',...
                                                        'Xtr_i');

    X_dot = reshape(Y(Xtr_i,:)',N_tr*n,1); 

    data = struct('X',X,'U',U,'X_dot',X_dot);
    
    load_PVTOL_params;
    
    %% Create Perm
    
    %do Perm matrix for training points
    Perm_Tr = zeros(N_tr*n,N*n);

    for i = 1:N_tr    
        idx = Xtr_i(i);
        Perm_Tr(1+(i-1)*n:i*n, 1+(idx-1)*n:idx*n) = eye(n);
    end
    Perm_Tr = sparse(Perm_Tr);

    features = struct('P_f',P_f,'P_b',P_b,'Perm_Tr',Perm_Tr);
    
    [mu_f_list, val_reg_err] = tune_L2_reg(X_val,U_val,Y_val,om_f,om_b,kernel_f.L,kernel_b.L,features,data,Xtr_i,N_tr,D_dyn);
    
    val_reg_err_all(:,nn) = val_reg_err;
    
end

%% Plot

close all;

linestyles = {'-*','-^','-s','-d'};
figure()
for nn = 1:4
    semilogy(mu_f_list,val_reg_err_all(:,nn),linestyles{nn},'linewidth',2); hold all
end
legend('N=100','N=250','N=500','N=1000');
xlabel('\mu_f'); ylabel('Mean reg. error');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

