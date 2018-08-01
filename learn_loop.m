
yalmip('clear');

%% Pre-compute features & jacobians

flag_sample_om = 1; %prevent accidental re-sample of features
fprintf('Computing features & Jacobians...');
[P_f,JP_f_all,JP_f,om_f] = Phi_f_T(X,O_dyn,Xc_i,kernel_f,m,n_dyn);
[P_b,om_b] = Phi_b_T(X,O_dyn,kernel_b,n_dyn_B);
[P_Wp_all,P_Wp, P_Wnp_all,P_Wnp,JP_Wp_all,JP_Wp,om_w_p,om_w_np] = Phi_w(n_W,X,Xc_i,O_w,kernel_wp,kernel_wnp,m,flag_sample_om,[],[]);
fprintf('Done!\n');
flag_sample_om = 0;

%% Start Loop

%termination conditions
max_itr = 15;
eps_term = 1e-1;

alpha_prev = zeros(D_dyn,1);
beta_prev = zeros(D_dyn,m);
lambda_prev = 0;
theta_p_prev = zeros(D_w,(n-m)*(n-m+1)/2);
theta_np_prev = zeros(D_w,n*(n+1)/2);

%Initialize max slack
max_slack = 100;

%Initialize change in solution
var_diff = 100;

%iteration count
itr = 0;

%keep track of constraint viol over entire set 
constraint_viol = nan(N,max_itr);

%% First do an independent unconstrained solve

%Unconstrainted, regularized baseline
[alpha_un, beta_un] = solve_dynamics_unconstrained(P_f,P_b,U,Y,mu_f,mu_b,n,m,Xtr_i,N_tr,D_dyn);
[f_h, B_h, df_h, dB_h] = construct_dyn(om_f, alpha_un, om_b, beta_un, O_dyn, kernel_f.L, kernel_b.L,n,m);
save(strcat('learned_functions/PVTOL_uncon_Dyn_Functions_',num2str(N_tr),'.mat'),...
        'f_h','B_h','df_h','dB_h');
clear f_h B_h df_h dB_h alpha_un beta_un;

%Unconstrained, un-regularized baseline
[alpha_unu, beta_unu] = solve_dynamics_unconstrained(P_f,P_b,U,Y,0,1e-6,n,m,Xtr_i,N_tr,D_dyn);
[f_h, B_h, df_h, dB_h] = construct_dyn(om_f, alpha_unu, om_b, beta_unu, O_dyn, kernel_f.L, kernel_b.L,n,m);
save(strcat('learned_functions/PVTOL_unconu_Dyn_Functions_',num2str(N_tr),'.mat'),...
        'f_h','B_h','df_h','dB_h');
clear f_h B_h df_h dB_h alpha_unu beta_unu;

%% Initialize metric

W = repmat(eye(n),1,1,length(Xc_i));
W_p = repmat(eye(n-m,n-m),1,1,length(Xc_i));
partial_Wp = zeros(n-m,n-m,length(Xc_i),n-m);

while (itr < max_itr) && (var_diff > eps_term)
    
    %% solve dynamics
    [alpha, beta, lambda] = solve_dynamics_opt(alpha_prev,beta_prev,P_f,JP_f,P_b,W_p,W,partial_Wp,U,Y,mu_f,mu_b,mu_s,n,m,N_tr,Xtr_i,Nc,Xc_i,D_dyn,eps_l,delta_l,max_slack);
    
    %% compute f(X), df(Xc)
    [f, dfdx_p] = compute_f_df(alpha,P_f,JP_f,n,m,Nc);
    
    %% solve W
    [theta_p, theta_np, w_lower] = solve_W_prob(theta_p_prev,theta_np_prev,P_Wp,P_Wnp,JP_Wp,n_W,lambda,f,dfdx_p,mu_w,n,m,N,Nc,Xc_i,D_w,eps_l,eps_wl,delta_wl,max_slack);
    
    %% update constraint set
    %compute df(X)
    [~, dfdx_p_all] = compute_f_df(alpha,P_f,JP_f_all,n,m,N);
    
    %compute W(X), dW(X)
    [W_p,W] = W_fnc(P_Wp_all*theta_p,P_Wnp_all*theta_np,n,m);
    partial_Wp = dWp_fnc(JP_Wp_all,theta_p,N,n,m);
    
    %update Xc
    [Xc_i,constraint_viol(:,itr+1),flag_done,max_slack] = update_constraint_set(eps_l,delta_wl,eps_wl,tol,W_p,W,partial_Wp,f,dfdx_p_all,lambda,N,Nc,n,m,Xc_i);
    Nc = length(Xc_i);
    
    %update relevant matrices
    if (flag_done) %all constraints satisfied.
        break;
    else
        JP_f = JP_f_all(:,:,Xc_i,:);
        P_Wp = P_Wp_all(Xc_i,:);
        P_Wnp = P_Wnp_all(Xc_i,:);
        JP_Wp = JP_Wp_all(Xc_i,:,:);
    end   
    
    %% update W(Xc), dW(Xc)
    [W_p,W] = W_fnc(P_Wp*theta_p,P_Wnp*theta_np,n,m);
    partial_Wp = dWp_fnc(JP_Wp,theta_p,Nc,n,m);
    
    %% check diff
    var_diff = max([norm(alpha-alpha_prev,'inf'),...
        norm(beta(:)-beta_prev(:),'inf'),...
        abs(lambda-lambda_prev),...
        norm(theta_p(:)-theta_p_prev(:),'inf'),...
        norm(theta_np(:)-theta_np_prev(:),'inf')]);
    
    itr = itr+1;
    fprintf('itr: %d, diff:%.4f \n',itr,var_diff);
    alpha_prev = alpha; beta_prev = beta; lambda_prev = lambda;
    theta_p_prev = theta_p; theta_np_prev = theta_np;
end

%% Define functions for f, B, df, W, dW

save(strcat('learned_functions/PVTOL_raw_',num2str(N_tr),'.mat'));

%NOT VECTORIZED
[f_h, B_h, df_h, dB_h] = construct_dyn(om_f, alpha, om_b, beta, O_dyn, kernel_f.L, kernel_b.L,n,m);
save(strcat('learned_functions/PVTOL_Dyn_Functions_',num2str(N_tr),'.mat'),...
        'f_h','df_h','B_h', 'dB_h');

%VECTORIZED
[W_h, dW_h] = construct_metric(om_w_p, theta_p, om_w_np, theta_np, O_w, n,m);
save(strcat('learned_functions/PVTOL_Metric_Functions_',num2str(N_tr),'.mat'),...
        'W_h','dW_h','lambda');

%% Evaluate regression quality & contraction constraints on full train & validation sets

run_validation(X,X_val,U_val,Y_val,alpha,beta,om_f,om_b,kernel_f.L,kernel_b.L,df_h,W_h,dW_h,lambda,eps_l,delta_wl,eps_wl);


