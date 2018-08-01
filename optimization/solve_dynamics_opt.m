function [alpha,beta,lambda] = solve_dynamics_opt(alpha_prev,beta_prev,Zf,Jf,Zb,W_p,W,partial_Wp,U,Y,mu_f,mu_b,mu_s,n,m,Ntr,Xtr_i,Nc,Xc_i,D,eps_l,delta_l,max_slack)

%Dynamics sub-problem

% yalmip('clear');

%% Define variables

%f:
alpha = sdpvar(D,1);

%b_j:
beta = sdpvar(D,m,'full');

%lambda:
lambda = sdpvar(1);

%error:
err = sdpvar(Ntr*n,1);

%slacks:
slacks = sdpvar(Nc,1);

%Define base constraints:
Constraints = [slacks >= 0;
               slacks <= max_slack*ones(Nc,1);
               lambda >= delta_l];

%% Define dynamics expressions
 
%f(X) (column-stacked)
f = Zf*alpha;

%B(X)*U
B = Zb*beta;

%restrict to training set
f_tr = repmat(alpha(1),Ntr*n,1);
Btr = repmat(beta(1),Ntr*n,m);

for i = 1:n
   f_tr(i:n:i+(Ntr-1)*n) = f(Xtr_i*n-(n-i)); 
end
for j = 1:m
    for i = 1:n
        Btr(i:n:i+(Ntr-1)*n,j) = B(Xtr_i*n-(n-i),j);
    end
end
Bu = Btr(:,1).*kron(U(Xtr_i,1),ones(n,1));
for j = 2:m
   Bu = Bu + Btr(:,j).*kron(U(Xtr_i,j),ones(n,1)); 
end

%reshape output
X_dot = reshape(Y(Xtr_i,:)',Ntr*n,1);

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
       dfdx_p(j,:) = alpha'*Jf(:,:,k,j); 
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
                                      +mu_s*sum(slacks)-lambda, options);    
 
fprintf('Done!');

%% Parse solution

if  (~sol.problem)
    alpha = double(alpha);
    beta = double(beta);
    lambda = double(lambda);
    slacks = double(slacks);
    err = reshape(double(err),Ntr,n);
    fprintf(' Regression err: %.3f, Max slack: %.3f. \n',(1/Ntr)*sum(norms(err,2,2)),max(slacks));
else
    fprintf('Problem!');
    keyboard;
end

end

    


