function [alpha,beta] = solve_dynamics_SDP(alpha_prev,beta_prev,features,data,constants,Xc_i,Nc,max_slack,...
                                                    J_Phi_f,W_p,W,partial_Wp)

%Dynamics sub-problem

yalmip('clear');

%% Get stuff

n = constants.n;
m = constants.m;
Ntr = constants.Ntr;
Xtr_i = constants.Xtr_i;
D = constants.D;

eps_l = constants.eps_l;

mu_f = constants.mu_f;
mu_b = constants.mu_b;
mu_s = constants.mu_s;

Phi_f = features.P_f;
Phi_b = features.P_b;
Perm_Tr = features.Perm_Tr;

U = data.U;
X_dot = data.X_dot;

%% Define variables

%f:
alpha = sdpvar(D,1);

%b_j:
beta = sdpvar(D,m,'full');

%lambda:
% lambda = sdpvar(1);
lambda = 0.1;

%error:
err = sdpvar(Ntr*n,1);

%slacks:
slacks = sdpvar(Nc,1);

%Define base constraints:
Constraints = [slacks >= 0;
               slacks <= max_slack*ones(Nc,1)];

%% Define dynamics expressions
 
f = Phi_f * alpha;

f_tr = Perm_Tr * f;
Bu = zeros(Ntr*n,1);
for j = 1:m
    Bu = Bu + sparse(diag(kron(U(Xtr_i,j),ones(n,1)))) * Perm_Tr * Phi_b * beta(:,j);
end

%Define error:
Constraints = [Constraints;
               err == f_tr + Bu - X_dot];

%% Define contraction constraints

%initialize
dfdx_p = repmat(alpha(1),n-m,n);

dWp_f_all = zeros(Nc*(n-m)*(n-m),1);
for j = 1:n-m
   f_j_Nc = repmat(f(Xc_i*n-(n-j)),1,(n-m)*(n-m));
   partial_Wp_j = partial_Wp(:,:,:,j);
   dWp_f_all = dWp_f_all + partial_Wp_j(:).*reshape(f_j_Nc',(n-m)*(n-m)*Nc,1);
end

for k = 1:Nc
    
    dWp_f = reshape(dWp_f_all(1+((n-m)*(n-m))*(k-1):k*((n-m)*(n-m))),(n-m),(n-m));
    
    for j = 1:n-m
       dfdx_p(j,:) = alpha'*J_Phi_f(:,:,k,j); 
    end
    
    Constraints = [Constraints;
            -dWp_f + dfdx_p*W(:,1:n-m,k) + W(1:n-m,:,k)*dfdx_p' <= -2*(lambda+eps_l)*W_p(:,:,k) + slacks(k)*eye(n-m)];
    
end

options = sdpsettings('solver','mosek','verbose',0);
% options = sdpsettings('solver','scs-direct','verbose',0);
% options.scs.max_iters=6000;
% options.scs.eps = 1.00e-06;
fprintf('Solving Dynamics problem...');
sol = optimize(Constraints, (err'*err)+...
                                      +mu_f*((alpha-alpha_prev)'*(alpha-alpha_prev))+...
                                      +mu_b*trace((beta-beta_prev)'*(beta-beta_prev))+...
                                      +mu_s*sum(slacks), options);    
 
fprintf('Done!');

%% Parse solution

if  (~sol.problem)
    alpha = double(alpha);
    beta = double(beta);
    slacks = double(slacks);
    err = reshape(double(err),n,Ntr);
    fprintf(' Regression err: %.3f, Max F viol: %.3f. \n',(1/Ntr)*sum(norms(err,2,1)),max(slacks));
else
    fprintf('Problem!');
    keyboard;
end

end

    


