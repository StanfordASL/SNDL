function [mu_f_list, val_reg_err] = tune_L2_reg(X_val,U_val,Y_val,om_f,om_b,Lf,Lb,features,data,Xtr_i,N_tr,D_dyn)


P_f = features.P_f;
P_b = features.P_b;
Perm_Tr = features.Perm_Tr;

U = data.U;
X_dot = data.X_dot;

%%

N_test = size(X_val,1);
n = size(X_val,2);
m = size(U_val,2);

O_f = size(om_f,1);

%% Define vectorized dynamics fncs for val dataset

Phi_T_f = (1/sqrt(O_f)) * ( kron ( kron(cos(X_val*om_f'),[1,0]) + kron(sin(X_val*om_f'),[0,1]) , Lf' ) );
Phi_T_b = (1/sqrt(O_f)) * ( kron ( kron(cos(X_val*om_b'),[1,0]) + kron(sin(X_val*om_b'),[0,1]) , Lb' ) );

X_dot_test = reshape(Y_val',N_test*n,1);


%% Reg constant tuning for Ridge-regression model

mu_b = 1e-6; %fixed
mu_f_list = [1e-3, 0.5e-3, 1e-4, 0.5e-4, 1e-5, 1e-6, 0.0];
val_reg_err = zeros(length(mu_f_list),1);

%Train & validate
for i = 1:length(mu_f_list)
    mu_f = mu_f_list(i);
%     mu_b = 10*mu_f;
    [alpha_un, beta_un] = solve_dynamics_unconstrained(Perm_Tr,P_f,P_b,U,X_dot,mu_f,mu_b,n,m,Xtr_i,N_tr,D_dyn);
    
    f_test = Phi_T_f*alpha_un;
    Bu_test = get_Bu(Phi_T_b,beta_un,U_val,n,m);
    
    reg_err = reshape(f_test + Bu_test - X_dot_test,n,N_test);
    val_reg_err(i) = (1/N_test)*sum(norms(reg_err,2,1)); 
end

end

function Bu = get_Bu(Phi_T_b,beta,U,n,m)
B_test = Phi_T_b*beta;
Bu = B_test(:,1).*kron(U(:,1),ones(n,1));
for j = 2:m
    Bu = Bu + B_test(:,j).*kron(U(:,j),ones(n,1));
end
end