function [Theta_p,Theta_np] = solve_metric(theta_p_prev,theta_np_prev,lambda,data,f,dfdX_p,constants,Xc_i,Nc,...
                                                    alpha,Zw_p,Zw_np,Jw_p,dF_dtheta_h,dW_h)


%% Get stuff

n = constants.n;
m = constants.m;
D = constants.D_w;
eps_l = constants.eps_l;
mu_s = constants.mu_s;
mu_w = constants.mu_w;
n_W = constants.n_W;
delta_wl = constants.delta_wl;
w_const = constants.w_const;
eps_wl = constants.eps_wl;
mu_w_phase1 = constants.mu_w_phase1;

optim = constants.optim;

X = data.X;

%% Problem variables

n_p = (n-m)*(n-m+1)/2;
n_np = (n*(n+1)/2) - n_p;

%theta_p : D x n_p
%theta_np: D x n_np

n_vars = D*(n_p + n_np);

%% Compute some problem info

dW_dtheta = zeros(n^2,n_vars,Nc);
dF_dtheta = zeros((n-m)^2,n_vars,Nc);

for k = 1:Nc
    
    x = X(Xc_i(k),:)';
    
    dW_dtheta_all = dW_h(x);
    dW_dtheta(:,:,k) = [dW_dtheta_all{1}, dW_dtheta_all{2}];
    
    dF_dtheta_all = dF_dtheta_h(x, alpha, lambda);
    dF_dtheta(:,:,k) = [dF_dtheta_all{1}, dF_dtheta_all{2}];
    
end

        
%% Create objective & grad for metric infeasibility

%for minimizing infeas
obj = @(vars,prob) metric_infeas_obj(vars,prob,theta_p_prev,theta_np_prev,Zw_p,Zw_np,Jw_p,f,dfdX_p,Xc_i,D,n,m,Nc,n_p,w_const,delta_wl,eps_wl,lambda,eps_l,constants.pen_tau);
JH = @(vars,prob) metric_infeas_grad(vars,prob,theta_p_prev,theta_np_prev,Nc,dW_dtheta,dF_dtheta,n,m,n_vars,constants.pen_tau);

%% Initial guess

var_guess = [theta_p_prev;
             theta_np_prev];
         
var_prev = var_guess;         

%% Minimize Metric uniform definiteness infeasibility

N_steps = optim.N;
term_tol = optim.term_tol;

%backtracking search params
line_alpha = optim.alpha;
line_beta = optim.beta;
t_lb = optim.t_lb;

prob_info = struct('mu_w_phase1',mu_w_phase1,'check_W',1);

eps_H = optim.eps_H*eye(n_vars);
phase_1 = 1;

while (phase_1)
    
    [vars,prob] = Newton_descent(var_guess,N_steps,term_tol,line_alpha,line_beta,t_lb,obj,JH,eps_H,prob_info);
    if (max(prob.MAX_W) < 0)
       phase_1 = 0;
    else
       mu_w_phase1 = 0.1*mu_w_phase1; 
       prob_info.mu_w_phase1 = mu_w_phase1;
       var_guess = vars;
    end
end

%% Setup const for metric problem

max_slack = max(prob.MAX_F)+0.1;

%% Solve problem

[theta,w_lower,w_upper,max_F,time] = solve_metric_SDP_single(var_prev,n_vars,D,n_p,Nc,n,m,n_W,delta_wl,eps_wl,w_const,lambda,eps_l,mu_w,mu_s,...
                                            Xc_i,Zw_p,Zw_np,Jw_p,f,dfdX_p,max_slack);


fprintf('w_lower: %.4f, w_upper: %.4f, max F slack: %.4f, solve time: %.2f \n', w_lower,w_upper,max_F,time);

Theta_p = reshape(theta(1:D*n_p),D,n_p);
Theta_np = reshape(theta(D*n_p+1:end),D,n_np);
    


end

    


