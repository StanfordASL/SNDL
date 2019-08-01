
flag_sample_om = 1; %prevent accidental re-sample of features

%% Pre-compute features & jacobians

fprintf('Computing features & Jacobians...');
[P_f,JP_f_all,JP_f,om_f] = Phi_f_T(X,O_dyn,Xc_i,kernel_f,m,n_dyn);
[P_b,om_b] = Phi_b_T(X,O_dyn,kernel_b,n_dyn_B);
[P_Wp_all,P_Wp, P_Wnp_all,P_Wnp,JP_Wp_all,JP_Wp,om_w_p,om_w_np] = Phi_w(n_W,X,Xc_i,O_w,kernel_wp,kernel_wnp,m,flag_sample_om,[],[]);
fprintf('Done!\n');
flag_sample_om = 0;


%do Perm matrix for training points
Perm_Tr = zeros(N_tr*n,N*n);

for i = 1:N_tr    
    idx = Xtr_i(i);
    Perm_Tr(1+(i-1)*n:i*n, 1+(idx-1)*n:idx*n) = eye(n);
end
Perm_Tr = sparse(Perm_Tr);

features = struct('P_f',P_f,'P_b',P_b,'Perm_Tr',Perm_Tr);
    
%% Get Function gradient handles (For step 2 Newton descent)

[dF_dalpha_h, dF_dtheta_h, dW_h] = gen_func_gradients(om_f,om_b,kernel_f,kernel_b,n,m,om_w_p,om_w_np,eps_l,w_const);

%% First do an independent unconstrained solve

%Unconstrainted, regularized baseline
[alpha_un, beta_un] = solve_dynamics_unconstrained(Perm_Tr,P_f,P_b,U,X_dot,1.0e-4,mu_b,n,m,Xtr_i,N_tr,D_dyn);
[f_h, B_h, df_h, dB_h] = construct_dyn(om_f, alpha_un, om_b, beta_un, O_dyn, kernel_f.L, kernel_b.L,n,m);
save(strcat('learned_functions/PVTOL_H_uncon_Dyn_Functions_',num2str(N_tr),'.mat'),...
        'f_h','B_h','df_h','dB_h');
clear f_h B_h df_h dB_h;

% %Unconstrainted, unregularized baseline
[alpha_unu, beta_unu] = solve_dynamics_unconstrained(Perm_Tr,P_f,P_b,U,X_dot,0.0,mu_b,n,m,Xtr_i,N_tr,D_dyn);
[f_h, B_h, df_h, dB_h] = construct_dyn(om_f, alpha_unu, om_b, beta_unu, O_dyn, kernel_f.L, kernel_b.L,n,m);
save(strcat('learned_functions/PVTOL_H_unconu_Dyn_Functions_',num2str(N_tr),'.mat'),...
        'f_h','B_h','df_h','dB_h');
clear f_h B_h df_h dB_h;

%% Start Loop

alpha_prev = zeros(D_dyn,1);
beta_prev = zeros(D_dyn,m);

n_p = (n-m)*(n-m+1)/2;
n_np = (n*(n+1)/2) - n_p;
Theta_p_prev = zeros(D_w,n_p);
Theta_np_prev = zeros(D_w,n_np);

%Initialize max slack
max_slack = 100;

%Initialize change in solution
var_diff = 100;

%iteration count
itr = 0;

%keep track of constraint viol over entire set 
constraint_viol = nan(N,max_itr);

%keep track of solutions
alpha_all = nan(D_dyn,max_itr);
beta_all = nan(D_dyn,m,max_itr);
theta_all = nan(D_w*(n_p+n_np),max_itr);

%% Initialize metric

[W_p,W] = W_fnc(P_Wp*Theta_p_prev,P_Wnp*Theta_np_prev,n,m,1); %initialize W as identity
partial_Wp = dWp_fnc(JP_Wp,Theta_p_prev,Nc,n,m);

%flag for done
flag_done = 0;


while (itr < max_itr) && (var_diff > eps_term) && (~flag_done)
    
    %% solve dynamics
    
    fprintf('***************************************************\n');
    fprintf('DYNAMICS\n');

    [alpha,beta] = solve_dynamics_SDP(alpha_prev,beta_prev,features,data,constants,Xc_i,Nc,max_slack,...
                                                    JP_f,W_p,W,partial_Wp);
                                                
    %% compute f(X), df(Xc)
    [f, dfdx_p] = compute_f_df(alpha,P_f,JP_f,n,m,Nc);
    
    %% solve W
    
    fprintf('***************************************************\n');
    fprintf('METRIC FEASIBILITY\n');
    [Theta_p, Theta_np] = solve_metric(reshape(Theta_p_prev,[],1),reshape(Theta_np_prev,[],1),lambda,data,f,dfdx_p,constants,Xc_i,Nc,...
                                                    alpha,P_Wp,P_Wnp,JP_Wp,dF_dtheta_h,dW_h);
    
    %% update constraint set
    %compute df(X)
    [~, dfdx_p_all] = compute_f_df(alpha,P_f,JP_f_all,n,m,N);
    
    %compute W(X), dW(X)
    [W_p,W] = W_fnc(P_Wp_all*Theta_p,P_Wnp_all*Theta_np,n,m,w_const);
    partial_Wp = dWp_fnc(JP_Wp_all,Theta_p,N,n,m);
    
    %update Xc
    fprintf('--------------------------------------------------\n');
    [Xc_i,constraint_viol(:,itr+1),flag_done,max_slack] = update_constraint_set(constants,W_p,W,partial_Wp,f,dfdx_p_all,lambda,Xc_i);
    Nc = length(Xc_i);
    
    %Update stuff at the new constraint points
    JP_f = JP_f_all(:,:,Xc_i,:);
    P_Wp = P_Wp_all(Xc_i,:);
    P_Wnp = P_Wnp_all(Xc_i,:);
    JP_Wp = JP_Wp_all(Xc_i,:,:);

    W_p = W_p(:,:,Xc_i);
    W = W(:,:,Xc_i);
    partial_Wp = partial_Wp(:,:,Xc_i,:);
    
    
    %% check diff
    var_diff = max([norm(alpha-alpha_prev,'inf'),...
        norm(beta(:)-beta_prev(:),'inf'),...
        norm(Theta_p(:)-Theta_p_prev(:),'inf'),...
        norm(Theta_np(:)-Theta_np_prev(:),'inf')]);
    
    itr = itr+1;
    fprintf('itr: %d, diff:%.4f \n',itr,var_diff);
    
    %% store solution
    
    alpha_all(:,itr) = alpha;
    beta_all(:,:,itr) = beta;
    theta_all(:,itr) = [Theta_p(:); Theta_np(:)];
    
    %% reset
    
    alpha_prev = alpha; beta_prev = beta;
    Theta_p_prev = Theta_p; Theta_np_prev = Theta_np;
    
end

%% Save

save(strcat('learned_functions/PVTOL_raw_',num2str(N_tr),'_new','.mat'));

%% Gen pareto

max_itr = itr;
reg_err_all = zeros(max_itr,1);
viol_all = zeros(max_itr,1);

for itr = 1:max_itr
    % Define functions
    [~, ~, df_h, ~] = construct_dyn(om_f, alpha_all(:,itr), om_b, beta_all(:,:,itr), O_dyn, kernel_f.L, kernel_b.L,n,m);
    [W_h, dWp_h] = construct_metric(om_w_p, reshape(theta_all(1:D_w*n_p,itr),D_w,n_p), ...
                                    om_w_np,reshape(theta_all(D_w*n_p+1:end,itr),D_w,n_np), O_w, n,m, w_const);    

    % Run validation test
    [~,~,~,...
          max_eig_F_test,~,min_eig_W_test,reg_err_val] = run_validation(X,X_val,U_val,Y_val,alpha_all(:,itr),beta_all(:,:,itr),...
                                        om_f,om_b,kernel_f.L,kernel_b.L,df_h,W_h,dWp_h,lambda,eps_l,delta_wl,eps_wl);

    reg_err_all(itr) = reg_err_val;
    
    constraint_viol_test = max([max_eig_F_test,(delta_wl+eps_wl)-min_eig_W_test],[],2);
    
    viol_all(itr) = mean(constraint_viol_test>0.0);
end

%% Plot pareto

% close all
figure()
legend_list = cell(max_itr,1);
hold all
for itr = 1:max_itr
    scatter(reg_err_all(itr),viol_all(itr),50,'filled'); 
    text(reg_err_all(itr),viol_all(itr)+0.01,sprintf('%d',itr));
    legend_list{itr} = sprintf('%d',itr);
end
plot(reg_err_all,viol_all,'ks-','linewidth',1.5,'markersize',15);
grid on; xlabel('Reg err Val'); ylabel('Fraction violations val');

keyboard;
%% Choose best itr and generate save functions


itr_best = max_itr;

alpha = alpha_all(:,itr_best);
beta = beta_all(:,:,itr_best);
Theta_p = reshape(theta_all(1:D_w*n_p,itr_best),D_w,n_p);
Theta_np = reshape(theta_all(D_w*n_p+1:end,itr_best),D_w,n_np);

%NOT VECTORIZED
[f_h, B_h, df_h, dB_h] = construct_dyn(om_f, alpha, om_b, beta, O_dyn, kernel_f.L, kernel_b.L,n,m);
save(strcat('learned_functions/PVTOL_Dyn_Functions_',num2str(N_tr),'_new','.mat'),...
        'f_h','df_h','B_h', 'dB_h');

%VECTORIZED
[W_h, dWp_h] = construct_metric(om_w_p,Theta_p,om_w_np,Theta_np, O_w, n,m, w_const); 
save(strcat('learned_functions/PVTOL_Metric_Functions_',num2str(N_tr),'_new','.mat'),...
        'W_h','dW_h','lambda');

%% Evaluate regression quality & contraction constraints on full train & validation sets

[max_eig_F_train,max_eig_W_train,min_eig_W_train,...
          max_eig_F_test,max_eig_W_test,min_eig_W_test,~] = run_validation(X,X_val,U_val,Y_val,alpha,beta,...
                                        om_f,om_b,kernel_f.L,kernel_b.L,df_h,W_h,dWp_h,lambda,eps_l,delta_wl,eps_wl);

%% Plots

close all

figure()
boxplot(constraint_viol);
title('Constraint viol vs iter');

figure()
hist_1 = histogram(max_eig_F_train,'facecolor','b','facealpha',.5,'Normalization','probability'); hold on
hist_2 = histogram(max_eig_F_test,'facecolor','r','facealpha',.5,'Normalization','probability');
hist_2.BinWidth = hist_1.BinWidth;
title('F violations');

figure()
subplot(1,2,1)
hist_1 = histogram(min_eig_W_train,'facecolor','b','facealpha',.5,'Normalization','probability'); hold on
hist_2 = histogram(min_eig_W_test,'facecolor','r','facealpha',.5,'Normalization','probability');
hist_2.BinWidth = hist_1.BinWidth;
title('min W spread');

subplot(1,2,2)
hist_1 = histogram(max_eig_W_train,'facecolor','b','facealpha',.5,'Normalization','probability'); hold on
hist_2 = histogram(max_eig_W_test,'facecolor','r','facealpha',.5,'Normalization','probability');
hist_2.BinWidth = hist_1.BinWidth;
title('max W spread');
