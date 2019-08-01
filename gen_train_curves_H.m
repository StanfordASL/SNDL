clear; close all; clc;

%% constraint_viol & 

CON_VIOL = nan(3100,15,2);
MAX_ITR = zeros(2,1);

REG_VIOL = cell(2,1);

%% load models

load('learned_functions/PVTOL_H_raw_150.mat');
max_itr = itr;
MAX_ITR(1) = max_itr;
CON_VIOL(:,1:max_itr,1) = constraint_viol(:,1:max_itr);
reg_err_all = zeros(itr,1);
viol_all = zeros(itr,1);
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
REG_VIOL{1} = [reg_err_all,viol_all];

clearvars -except CON_VIOL MAX_ITR REG_VIOL

%%

load('learned_functions/PVTOL_H_raw_1000_new.mat'); 
max_itr = itr;
MAX_ITR(2) = max_itr;
CON_VIOL(:,1:max_itr,2) = constraint_viol(:,1:max_itr);
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
REG_VIOL{2} = [reg_err_all,viol_all];

clearvars -except CON_VIOL MAX_ITR REG_VIOL

%% compute sigma, mean
mean_vec = nan(2,15);
sigma_vec = nan(2,15);

for i = 1:2
    mean_vec(i,1:MAX_ITR(i)) = mean(CON_VIOL(:,1:MAX_ITR(i),i));
    sigma_vec(i,1:MAX_ITR(i)) = std(CON_VIOL(:,1:MAX_ITR(i),i));
end

%% cbrewer

cmap = cbrewer('qual','Set1',2);

%% gen plot

[~,order] = sort(mean_vec(:,1)+sigma_vec(:,1),'descend');
linestyles = {'-*','-^','-s','-d'};
close all
figure()
hold on
for idx = 1:2
    i = order(idx);

    x = (1:MAX_ITR(i))';
    y = mean_vec(i,1:MAX_ITR(i))';
    dy = 2.5*sigma_vec(i,1:MAX_ITR(i))';
    
    fill([x;flipud(x)],[y-dy;flipud(y+dy)],cmap(i,:),'linestyle','none','facealpha',0.8);
    xlim([1,max(MAX_ITR)])
end

grid on
xlabel('Global iteration');
ylabel('\pm 2.5\sigma spread of \nu^{(k)}');
nom_legend = {'N = 150','N = 1000'};
for idx = 1:2
    nom_legend{idx} = sprintf('%s ; %s',nom_legend{idx},...
                                        strcat('#itr = ',num2str(MAX_ITR(idx))));
end
legend(nom_legend{order});
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

keyboard;
%% reg plot

close all;
figure()
hold all
for idx = 1:2
    plot(REG_VIOL{idx}(:,1),REG_VIOL{idx}(:,2),linestyles{idx},'linewidth',2,'markersize',12);
end
grid on
xlabel('Mean reg. error');
ylabel('Frac. violations');
legend('N = 150','N = 1000')
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

                             