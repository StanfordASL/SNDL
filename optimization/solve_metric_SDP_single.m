function [theta,w_lower,w_upper,max_F,time] = solve_metric_SDP_single(theta_prev,n_vars,D,n_p,Nc,n,m,n_W,delta_wl,eps_wl,w_const,lambda,eps_l,mu_w,mu_s,...
                                                Xc_i,Zw_p,Zw_np,Jw_p,f,dfdX_p,max_slack)

%metric problem

% yalmip('clear');

%% Define variables

theta = sdpvar(n_vars,1);

theta_p = theta(1:D*n_p);
theta_np = theta(D*n_p+1:end);

Theta_p = reshape(theta_p,D,[]);
Theta_np = reshape(theta_np,D,[]);

%slacks:
slacks = sdpvar(Nc,1);

%lower, upper
w_lower = sdpvar(1);
w_upper = sdpvar(1);

%Define base constraints:
Constraints = [w_lower >= delta_wl;
               slacks >= 0;
               slacks <= max_slack*ones(Nc,1)];

%% Compute W_p, W, partial_Wp

[W_p,W] = W_fnc(Zw_p*Theta_p,Zw_np*Theta_np,n,m,w_const);
partial_Wp = dWp_fnc(Jw_p,Theta_p,Nc,n,m);

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
            -dWp_f + dfdX_p(:,:,k)*W(:,1:n-m,k) + W(1:n-m,:,k)*dfdX_p(:,:,k)' <= -2*(lambda+eps_l)*W_p(:,:,k) + slacks(k)*eye(n-m)];

                     
    Constraints_PD = [Constraints_PD;
                     (w_lower+eps_wl)*eye(n) <= W(:,:,k);
                      W(:,:,k) <= w_upper*eye(n)];
    
end

Constraints = [Constraints;
               Constraints_F;
               Constraints_PD];

%%

options = sdpsettings('solver','mosek','verbose',1);
% options = sdpsettings('solver','scs-direct','verbose',1);
% options.scs.max_iters=10000;
% options.scs.eps = 1.00e-06;
% fprintf('Solving metric problem...');
sol = optimize(Constraints, (mu_s/Nc)*(w_upper-w_lower)+(mu_s/Nc)*mu_w*((theta-theta_prev)'*(theta-theta_prev))+...
                                             +(1/Nc)*sum(slacks),options); 
% fprintf('Done!');

%% Parse solution

if  (~sol.problem)

    theta = double(theta);   
    w_lower = double(w_lower);
    w_upper = double(w_upper);
    slacks = double(slacks);
    max_F = max(slacks);
    time = sol.solvertime;
%     fprintf(' w_lower = %.3f, w_upper = %.3f, Max slack: %.3f \n',w_lower, w_upper, max(slacks));
else
    fprintf('Problem!');
    keyboard;
end

end
