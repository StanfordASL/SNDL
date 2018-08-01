function [alpha,beta] = solve_dynamics_unconstrained(Zf,Zb,U,Y,mu_f,mu_b,n,m,Xtr_i,Ntr,D)

%unconstrained dynamics fitting (with or without l2 reg)

%% Define variables

%f:
alpha = sdpvar(D,1);

%b_j:
beta = sdpvar(D,m,'full');

%error:
err = sdpvar(Ntr*n,1);


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
Constraints = [err == f_tr + Bu - X_dot];

options = sdpsettings('solver','mosek','verbose',0);
% options = sdpsettings('solver','scs-direct','verbose',0);
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
    err = reshape(double(err),Ntr,n);
    fprintf('Regression error : %.3f\n',(1/Ntr)*sum(norms(err,2,2)));
end

end

    


