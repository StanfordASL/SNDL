function [dF_dalpha_h, dF_dtheta_h, dW_h] = gen_func_gradients(om_f,om_b,kernel_f,kernel_b,n,m,om_w_p,om_w_np,eps_l,w_const)
% function dF_dalpha_h = gen_func_gradients(N_tr)
    
%% functional dependence

% f_h (x, alpha)
% B_h (x, beta)
 
% W_h (x, theta_p, theta_np)
 
% F_h (x,alpha,theta_p,theta_np,lambda)
% dF_h {alpha, theta_p, theta_np, lambda}


%% Import random features

% load(sprintf('learned_functions/PVTOL_raw_%d',N_tr), 'om_f','om_b','O_dyn','kernel_f','kernel_b','n','m','om_w_p','om_w_np','eps_l',...
%                                                       'mu_f','mu_b',);

%%

% Identify some constants
O_f = size(om_f, 1);        % number of Gaussian directions for f(x)
O_b = size(om_b, 1);        % ... B(x)
D_f = 2 * O_f * n;          % dimension of parameter vector for f(x)
D_b = 2 * O_b * n;          % ... B(x)

O_w_p = size(om_w_p, 1);
O_w_np = size(om_w_np, 1);
D_w_p = 2 * O_w_p;
D_w_np = 2 * O_w_np;

Lf = kernel_f.L;
Lb = kernel_b.L;

%% Instantiate symbolic variables

x = tom('x', n, 1);
alpha = tom('alpha', D_f, 1);
beta = tom('beta', D_b, m);

lambda = tom('lambda');

num_p = (n-m)*(n-m+1)/2;
num_np = n*(n+1)/2 - num_p;

theta_p = tom('theta_p',D_w_p*num_p,1);
theta_np = tom('theta_np', D_w_np*num_np,1);

Theta_p = reshape(theta_p,D_w_p,num_p);
Theta_np = reshape(theta_np, D_w_np,num_np);

%% Kernel features
phi_f_T = (1/sqrt(O_f))*(kron(kron(cos(x'*om_f'), [1,0]) + kron(sin(x'*om_f'), [0,1]), Lf'));               % n x D_f
phi_b_T = (1/sqrt(O_b))*(kron(kron(cos(x'*om_b'), [1,0]) + kron(sin(x'*om_b'), [0,1]), Lb'));               % n x D_b
phi_w_p_T = (1/sqrt(O_w_p))*(kron(cos(x(1:n-m)'*om_w_p'), [1,0]) + kron(sin(x(1:n-m)'*om_w_p'), [0,1]));    % 1 x D_w_p
phi_w_np_T = (1/sqrt(O_w_np))*(kron(cos(x'*om_w_np'), [1,0]) + kron(sin(x'*om_w_np'), [0,1]));              % 1 x D_w_np

%% Dynamics functions (symbolic)
f = phi_f_T * alpha;                            % f(x, alpha),              n x 1
B = phi_b_T * beta;                             % B(x, beta),               n x m


% dfdx_p
Jphi_j_f = (1/sqrt(O_f))*(kron(-diag(sin(om_f*x))*om_f,[1;0]) + kron(diag(cos(om_f*x))*om_f,[0;1]));
% dfdx_p = tom('dfdx_p',n-m,n);

for i = 1:n-m
    dfdx_p(i,:) = alpha'*kron(Jphi_j_f,Lf(:,i));
end
           
%% Metric function (symbolic)

W_l = w_const*tomSym(eye(n));

W_rest = vec_to_sym_tom([phi_w_p_T*Theta_p, phi_w_np_T*Theta_np],n);

W =  W_l + W_rest;

W_p = W(1:(n-m),1:(n-m));

%% Contraction LMI

%dphi/dx_j (X):
Jw_p_j = @(j) (1/sqrt(O_w_p))* (kron(-sin(x(1:n-m)'*om_w_p').*om_w_p(:,j)',[1,0])+...
                                kron( cos(x(1:n-m)'*om_w_p').*om_w_p(:,j)',[0,1]));
                              
dWp_f = f(1)*Jw_p_j(1)*Theta_p;
for j = 2:n-m
    dWp_f = dWp_f + f(j)*Jw_p_j(j)*Theta_p;
end

dWp_f = vec_to_sym_tom(dWp_f,n-m);

F = -dWp_f + dfdx_p*W(:,1:n-m) + W(1:n-m,:)*dfdx_p' + 2*(lambda+eps_l)*W_p;

%% Gradients

dF_dalpha = derivative(F,alpha);

dF_dtheta_p = derivative(F,theta_p);

dF_dtheta_np = derivative(F,theta_np);

dW_dtheta_p = derivative(W,theta_p);

dW_dtheta_np = derivative(W, theta_np);

%% Convert to matlab function

dir_name = './function_gen_eval/tomTemp_PVTOL';
% mkdir(dir_name);

% dF_dalpha_h(x, theta_p, theta_np)
tempD = mfile(dF_dalpha, 'dF_dalpha_h_tom', dir_name);
dF_dalpha_h = @(x, theta_p, theta_np) dF_dalpha_h_tom(theta_np, theta_p, x, tempD);

% dF_dtheta_p_h(x, alpha, lambda)
tempD = mfile(dF_dtheta_p, 'dF_dtheta_p_h_tom', dir_name);
dF_dtheta_p_h = @(x, alpha, lambda) dF_dtheta_p_h_tom(alpha,lambda,x,tempD);

% dF_dtheta_np_h(x, alpha, lambda)
tempD = mfile(dF_dtheta_np, 'dF_dtheta_np_h_tom', dir_name);
dF_dtheta_np_h = @(x, alpha, lambda) dF_dtheta_np_h_tom(alpha,lambda,x,tempD);

dF_dtheta_h = @(x, alpha, lambda) { dF_dtheta_p_h(x, alpha, lambda) , dF_dtheta_np_h(x, alpha, lambda) };

% dW_dtheta_p(x)
tempD = mfile(dW_dtheta_p, 'dW_dtheta_p_h_tom', dir_name);
dW_dtheta_p_h = @(x) dW_dtheta_p_h_tom(x,tempD);

% dW_dtheta_np(x)
tempD = mfile(dW_dtheta_np, 'dW_dtheta_np_h_tom', dir_name);
dW_dtheta_np_h = @(x) dW_dtheta_np_h_tom(x,tempD);


dW_h = @(x) {dW_dtheta_p_h(x), dW_dtheta_np_h(x) };

end

function V = vec_to_sym_tom(v,n)

%symmetric matrix from minimal encoding vector
%V: n x n x N stack of symmetric matrices
%v: N x n(n+1)/2 stack of row vectors

%pre-allocate matrix
% V = repmat(v(1),n,n);

%reshape v
v = reshape(v', 1,size(v,2));
vp = v'; %permute(v,[2,1]);

%top-left corner
V(1,1) = v(1,1);

%indices
i_end = cumsum(1:n);

for j = 2:n
    V(j,1:j) = v(1,i_end(j-1)+1:i_end(j)); %row
end

for j = 2:n
    V(1:j-1,j) = vp(i_end(j-1)+1:i_end(j)-1,1); %column
end


end

function V = vec_to_tri_tom(v,n)

%symmetric matrix from minimal encoding vector
%V: n x n x N stack of symmetric matrices
%v: N x n(n+1)/2 stack of row vectors

%pre-allocate matrix
% V = repmat(v(1),n,n);

%reshape v
v = reshape(v', 1,size(v,2));
vp = v'; %permute(v,[2,1]);

%top-left corner
V(1,1) = v(1,1);

%indices
i_end = cumsum(1:n);

for j = 2:n
    V(j,1:j) = v(1,i_end(j-1)+1:i_end(j)); %row
end


end
