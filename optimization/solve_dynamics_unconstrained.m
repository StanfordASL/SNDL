function [alpha,beta] = solve_dynamics_unconstrained(Perm_Tr,Zf,Zb,U,X_dot,mu_f,mu_b,n,m,Xtr_i,Ntr,D)

%unconstrained dynamics fitting (with or without l2 reg)

%% Define variables

%f:
alpha = sdpvar(D,1);

%b_j:
beta = sdpvar(D,m,'full');

%error:
err = sdpvar(Ntr*n,1);


%% Define dynamics expressions
 
f_tr = Perm_Tr * Zf * alpha;
Bu = zeros(Ntr*n,1);
for j = 1:m
    Bu = Bu + diag(kron(U(Xtr_i,j),ones(n,1))) * Perm_Tr * Zb * beta(:,j);
end

%Define error:
Constraints = [err == f_tr + Bu - X_dot];

options = sdpsettings('solver','mosek','verbose',0);
% options = sdpsettings('solver','scs-direct','verbose',1);
% options.scs.max_iters=6000;
% options.scs.eps = 1.00e-06;
% tic
fprintf('Solving unconstrained dynamics problem...');
soln_unconstrained = optimize(Constraints, (err'*err)+...
                               +mu_f*(alpha'*alpha)+...
                               +mu_b*trace(beta'*beta), options); 
% elapsed=toc;
fprintf('Done!');

if soln_unconstrained.problem == 0
    alpha = double(alpha);
    beta = double(beta);
    err = reshape(double(err),n,Ntr);
    fprintf('Regression error : %.3f\n',(1/Ntr)*sum(norms(err,2,1)));
end

end

    


