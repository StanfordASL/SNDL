function [theta_p,theta_np,w_lower] = solve_W_prob(theta_p_prev,theta_np_prev,Zw_p,Zw_np,Jw_p,n_W,lambda,f,dfdz_p,mu_w,n,m,N,Nc,Xc_i,D,eps_l,eps_wl,delta_wl,max_slack)

%metric problem

% yalmip('clear');

%% Define variables

%Wp
theta_p = sdpvar(D,(n-m)*(n-m+1)/2,'full');

%Wnp
theta_np = sdpvar(D,n*(n+1)/2,'full');

%slacks:
slacks = sdpvar(Nc,1);

%lower, upper
w_lower = sdpvar(1);
w_upper = sdpvar(1);

%Define base constraints:
Constraints = [theta_np(:,1:(n-m)*(n-m+1)/2) == 0; %top (n-m)x(n-m) block zero
               w_lower >= delta_wl;
               slacks >= 0;
               slacks <= max_slack*ones(Nc,1)];

%% Compute W_p, W, partial_Wp

[W_p,W] = W_fnc(Zw_p*theta_p,Zw_np*theta_np,n,m);
partial_Wp = dWp_fnc(Jw_p,theta_p,Nc,n,m);

%% Define contraction & uniform definiteness constraints

dWp_f_all = zeros(Nc*(n-m)*(n-m),1);
for j = 1:n-m
   if ismember(j,n_W)
       continue;
   end
   f_j_Nc = repmat(f(Xc_i*n-(n-j)),1,(n-m)*(n-m));
   partial_Wp_j = partial_Wp(:,:,:,j);
   dWp_f_all = dWp_f_all + partial_Wp_j(:).*reshape(f_j_Nc',(n-m)*(n-m)*Nc,1);
end

Constraints_F = [];
Constraints_PD = [];

for k = 1:Nc

    dWp_f = reshape(dWp_f_all(1+((n-m)*(n-m))*(k-1):k*((n-m)*(n-m))),(n-m),(n-m));
       
    Constraints_F = [Constraints_F; 
            -dWp_f + dfdz_p(:,:,k)*W(:,1:n-m,k) + W(1:n-m,:,k)*dfdz_p(:,:,k)' <= -2*(lambda+eps_l)*W_p(:,:,k) + slacks(k)*eye(n-m)];
        
    Constraints_PD = [Constraints_PD;
                     (w_lower+eps_wl)*eye(n) <= W(:,:,k);
                      W(:,:,k) <= w_upper*eye(n)];
    
end

Constraints = [Constraints;
               Constraints_F;
               Constraints_PD];

options = sdpsettings('solver','mosek','verbose',0);
% options = sdpsettings('solver','scs-direct','verbose',0);
% options.scs.max_iters=6000;
% options.scs.eps = 1.00e-06;
fprintf('Solving metric problem...');
sol = optimize(Constraints, (w_upper-w_lower)+...
                                             +mu_w*trace((theta_p-theta_p_prev)'*(theta_p-theta_p_prev))+...
                                             +mu_w*trace((theta_np-theta_np_prev)'*(theta_np-theta_np_prev))+...
                                             +100*sum(slacks),options); 
fprintf('Done!');

%% Parse solution

if  (~sol.problem)
    theta_p = double(theta_p);
    theta_np = double(theta_np);
    slacks = double(slacks);
    w_lower = double(w_lower);
    w_upper = double(w_upper);
    fprintf(' w_lower = %.3f, w_upper = %.3f, Max slack: %.3f \n',w_lower, w_upper, max(slacks));
else
    fprintf('Problem!');
    keyboard;
end

end
