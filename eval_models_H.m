clear; close all; clc;

%% eval new models

load('PVTOL_H_raw_150.mat')

%% CCM model init
max_itr = itr;
itr_best = max_itr;

alpha_new = alpha_all(:,itr_best);
beta_new = beta_all(:,:,itr_best);
Theta_p_new = reshape(theta_all(1:D_w*n_p,itr_best),D_w,n_p);
Theta_np_new = reshape(theta_all(D_w*n_p+1:end,itr_best),D_w,n_np);

%construct metric from CCM model
[W_h, dWp_h] = construct_metric(om_w_p,Theta_p_new,om_w_np,Theta_np_new, O_w, n,m, w_const); 

%% Eval RR model

disp('RR');
[alpha_un, beta_un] = solve_dynamics_unconstrained(Perm_Tr,P_f,P_b,U,X_dot,mu_f,mu_b,n,m,Xtr_i,N_tr,D_dyn);
[f_h_un, B_h_un, df_h_un, dB_h_un] = construct_dyn(om_f, alpha_un, om_b, beta_un, O_dyn, kernel_f.L, kernel_b.L,n,m);
[max_eig_F_train,max_eig_W_train,min_eig_W_train,...
          max_eig_F_test,max_eig_W_test,min_eig_W_test,~] = run_validation(X,X_val,U_val,Y_val,alpha_un,beta_un,...
                                        om_f,om_b,kernel_f.L,kernel_b.L,df_h_un,W_h,dWp_h,lambda,eps_l,delta_wl,eps_wl);

fprintf('******************\n');

%% Eval CCM-R model

[f_h, B_h, df_h, dB_h] = construct_dyn(om_f, alpha_new, om_b, beta_new, O_dyn, kernel_f.L, kernel_b.L,n,m);

disp('CCM-R-new');
f_tr = Perm_Tr * P_f * alpha_new;
Bu = zeros(N_tr*n,1);
for j = 1:m
Bu = Bu + diag(kron(U(Xtr_i,j),ones(n,1))) * Perm_Tr * P_b * beta_new(:,j);
end
err = f_tr + Bu - X_dot;

fprintf('CCM train reg err: %.4f\n', (1/N_tr)*sum(norms(err,2,1)));

[max_eig_F_train,max_eig_W_train,min_eig_W_train,...
          max_eig_F_test,max_eig_W_test,min_eig_W_test,~] = run_validation(X,X_val,U_val,Y_val,alpha_new,beta_new,...
                                        om_f,om_b,kernel_f.L,kernel_b.L,df_h,W_h,dWp_h,lambda,eps_l,delta_wl,eps_wl);

constraint_viol_test = max([max_eig_F_test,(delta_wl+eps_wl)-min_eig_W_test],[],2);
constraint_viol_train = max([max_eig_F_train,(delta_wl+eps_wl)-min_eig_W_train],[],2);


fprintf('frac viol train: %.4f, frac viol test: %.4f \n',mean(constraint_viol_train>0.0),...
                                                            mean(constraint_viol_test>0.0) );
              
fprintf('******************\n');
