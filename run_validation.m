function run_validation(X_train,X_test,U_test,Y_test,alpha,beta,om_f,om_b,Lf,Lb,df_h,W_h,dW_h,lambda,eps_l,delta_wl,eps_wl)

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
      
W = W_h(X_train);
dW = dW_h(X_train);

max_eig_F = zeros(N_train,1);
min_eig_W = zeros(N_train,1);
max_eig_W = zeros(N_train,1);

Phi_T_f_train = (1/sqrt(O_f)) * ( kron ( kron(cos(X_train*om_f'),[1,0]) + kron(sin(X_train*om_f'),[0,1]) , Lf' ) );
f_train = Phi_T_f_train*alpha;

% Evaluate contraction constraint over training set

for i = 1:N_train
    
    f_i = f_train(1+(i-1)*n:i*n);
    dWp_f = zeros(n-m,n-m);
    for j = 1:n-m
        dWp_f = dWp_f + dW(1:n-m,1:n-m,i,j)*f_i(j);
    end
    
    dfdx_p = B_perp'*df_h(X_train(i,:)');
    
    F = -dWp_f + dfdx_p*W(:,1:n-m,i) + W(1:n-m,:,i)*dfdx_p' + 2*(lambda+eps_l)*W(1:n-m,1:n-m,i);
    
    max_eig_F(i) = max(eig(F));  
    min_eig_W(i) = min(eig(W(:,:,i)));
    max_eig_W(i) = max(eig(W(:,:,i)));
    
end

fprintf('************************\n');
fprintf('train: max contraction violation: %.3f \n', max(max_eig_F));
fprintf('train: min(eig(W)): %.3f (%.3f), max(eig(W)): %.3f \n',min(min_eig_W), delta_wl+eps_wl, max(max_eig_W));

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

reg_err = reshape(f_test + Bu - X_dot_test,N_test,n);

fprintf('test: reg_err = %.3f \n', (1/N_test)*sum(norms(reg_err,2,2)));

%% Check constraction constraints & uniform definiteness over validation/test set

W = W_h(X_test);
dW = dW_h(X_test);

max_eig_F = zeros(N_test,1);
min_eig_W = zeros(N_test,1);
max_eig_W = zeros(N_test,1);

for i = 1:N_test
    
    f_i = f_test(1+(i-1)*n:i*n);
    dWp_f = zeros(n-m,n-m);
    for j = 1:n-m
        dWp_f = dWp_f + dW(1:n-m,1:n-m,i,j)*f_i(j);
    end
    
    dfdx_p = B_perp'*df_h(X_test(i,:)');
    
    F = -dWp_f + dfdx_p*W(:,1:n-m,i) + W(1:n-m,:,i)*dfdx_p' + 2*(lambda+eps_l)*W(1:n-m,1:n-m,i);
    
    max_eig_F(i) = max(eig(F));  
    min_eig_W(i) = min(eig(W(:,:,i)));
    max_eig_W(i) = max(eig(W(:,:,i)));
    
end

fprintf('test: max contraction violation: %.3f, avg violation: %.3f \n', max(max_eig_F),mean(max_eig_F));
fprintf('test: min(eig(W)): %.3f (%.3f), max(eig(W)): %.3f \n',min(min_eig_W), delta_wl+eps_wl, max(max_eig_W));


end
