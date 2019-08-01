function [max_eig_F_train,max_eig_W_train,min_eig_W_train,...
          max_eig_F_test,max_eig_W_test,min_eig_W_test,reg_err_val] = run_validation(X_train,X_test,U_test,Y_test,alpha,beta,om_f,om_b,Lf,Lb,df_h,W_h,dW_h,lambda,eps_l,delta_wl,eps_wl)

%import full training and test dataset

N_train = size(X_train,1);
N_test = size(X_test,1);

n = size(X_test,2);
m = size(U_test,2);

O_f = size(om_f,1);

%% Check contraction constraints & uniform definiteness over all training set

%F: -partial_f Wp + dfdx_p*W*B_perp + B_perp'*W*dfdx_p + 2*lambda*Wp: negative
%definite

B_perp = [eye(n-m);
          zeros(m,n-m)];
      
W_train = W_h(X_train);
dW_p_train = dW_h(X_train);

max_eig_F_train = zeros(N_train,1);
min_eig_W_train = zeros(N_train,1);
max_eig_W_train = zeros(N_train,1);

Phi_f_T_train = (1/sqrt(O_f)) * ( kron ( kron(cos(X_train*om_f'),[1,0]) + kron(sin(X_train*om_f'),[0,1]) , Lf' ) );
f_train = Phi_f_T_train*alpha;

% Evaluate contraction constraint over training set

for i = 1:N_train
    
    f_i = f_train(1+(i-1)*n:i*n);
    dWp_f = zeros(n-m,n-m);
    for j = 1:n-m
        dWp_f = dWp_f + dW_p_train(:,:,i,j)*f_i(j);
    end
    
    dfdx_p = B_perp'*df_h(X_train(i,:)');
    
    F = -dWp_f + dfdx_p*W_train(:,1:n-m,i) + W_train(1:n-m,:,i)*dfdx_p' + 2*(lambda+eps_l)*W_train(1:n-m,1:n-m,i);
    
    max_eig_F_train(i) = max(eig(F));  
    min_eig_W_train(i) = min(eig(W_train(:,:,i)));
    max_eig_W_train(i) = max(eig(W_train(:,:,i)));
    
end

fprintf('************************\n');
fprintf('train: max contraction violation: %.3f \n', max(max_eig_F_train));
fprintf('train: min(eig(W)): %.3f (%.3f), max(eig(W)): %.3f \n',min(min_eig_W_train), delta_wl+eps_wl, max(max_eig_W_train));

%% Define vectorized dynamics fncs for test dataset

Phi_T_f = (1/sqrt(O_f)) * ( kron ( kron(cos(X_test*om_f'),[1,0]) + kron(sin(X_test*om_f'),[0,1]) , Lf' ) );
f_test = Phi_T_f*alpha;

Phi_T_b = (1/sqrt(O_f)) * ( kron ( kron(cos(X_test*om_b'),[1,0]) + kron(sin(X_test*om_b'),[0,1]) , Lb' ) );
B_test = Phi_T_b*beta;
Bu = B_test(:,1).*kron(U_test(:,1),ones(n,1));
for j = 2:m
   Bu = Bu + B_test(:,j).*kron(U_test(:,j),ones(n,1)); 
end

X_dot_test = reshape(Y_test',N_test*n,1);

%% Regression error over test set

reg_err = reshape(f_test + Bu - X_dot_test,n,N_test);
reg_err_val = (1/N_test)*sum(norms(reg_err,2,1));

fprintf('test: reg_err = %.3f \n',reg_err_val);

%% Check constraction constraints & uniform definiteness over validation/test set

W_test = W_h(X_test);
dW_p_test = dW_h(X_test);

max_eig_F_test = zeros(N_test,1);
min_eig_W_test = zeros(N_test,1);
max_eig_W_test = zeros(N_test,1);

for i = 1:N_test
    
    f_i = f_test(1+(i-1)*n:i*n);
    dWp_f = zeros(n-m,n-m);
    for j = 1:n-m
        dWp_f = dWp_f + dW_p_test(:,:,i,j)*f_i(j);
    end
    
    dfdx_p = B_perp'*df_h(X_test(i,:)');
    
    F = -dWp_f + dfdx_p*W_test(:,1:n-m,i) + W_test(1:n-m,:,i)*dfdx_p' + 2*(lambda+eps_l)*W_test(1:n-m,1:n-m,i);
    
    max_eig_F_test(i) = max(eig(F));  
    min_eig_W_test(i) = min(eig(W_test(:,:,i)));
    max_eig_W_test(i) = max(eig(W_test(:,:,i)));
    
end

fprintf('test: max contraction violation: %.3f, avg violation: %.3f \n', max(max_eig_F_test),mean(max_eig_F_test));
fprintf('test: min(eig(W)): %.3f (%.3f), max(eig(W)): %.3f \n',min(min_eig_W_test), delta_wl+eps_wl, max(max_eig_W_test));


end
