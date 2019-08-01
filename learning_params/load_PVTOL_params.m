%% Features

O_dyn = 8*n; %# random gaussian directions for dynamics
O_w = 6*n; %# random gaussian directions for W

D_dyn = 2*O_dyn*n;
D_w = 2*O_w;

%exclude these components of x in f(x)
n_dyn = 1;

%exclude these components of x in B(x)
n_dyn_B = [1,3];

%exclude these components of x in W(x)
n_W = 1;

% Kernel definitions
kernel_f = struct('sigma',6,'L',eye(n));
kernel_b = struct('sigma',6.0,'L',blkdiag(zeros(n-m),eye(m)));
kernel_wp = struct('sigma',15.0);
kernel_wnp = struct('sigma',15.0);

%% Regularization & tuning constants

mu_f = 1e-3;
mu_b = 10*mu_f;
mu_w = 1e-3;

mu_s = 1e-2;

mu_w_phase1 = mu_w;

lambda = 0.1;

%% Tolerance params

% slack tolerance for re-sample
tol = 0.05;

% extra tolerance for LMI constraints
eps_l = 0.1; %lambda
eps_wl = 0.1; %wl

% lower bounds on LMI constraints
delta_wl = 0.1; %wl

% manual extra constant for W
w_const = 0.0;

%% Optimization params

%termination conditions
max_itr = 15; %global iterations
eps_term = 1e-2;

%penalty params (to define penalty fncs)
pen_scale = 1;
pen_tau = -0.1;

%params for step 2 Newton descent
N_steps = 100; %max iters
term_tol = 1e-3; %termination tolerance
line_alpha = 0.01; %backtrack alpha
line_beta = 0.5; %backtrack beta
t_lb = 1e-5; %Newton step lower-bound
eps_H = (1.0e-10); %Hessian damping

constants.optim = struct('N',N_steps,'term_tol',term_tol,...
                         'alpha',line_alpha,'beta',line_beta,'t_lb',t_lb,...
                         'eps_H',eps_H);

%% Structify

constants.eps_term = eps_term;
constants.tol = tol;
constants.n = n;
constants.m = m;
constants.N = N;
constants.Ntr = N_tr;
constants.Xtr_i = Xtr_i;
constants.D = D_dyn;
constants.D_w = D_w;
constants.eps_l = eps_l;
constants.mu_f = mu_f;
constants.mu_b = mu_b;
constants.mu_s = mu_s;
constants.mu_w = mu_w;
constants.mu_w_phase1 = mu_w_phase1;
constants.n_W = n_W;
constants.delta_wl = delta_wl;
constants.eps_wl = eps_wl;
constants.w_const = w_const;
constants.pen_scale = pen_scale;
constants.pen_tau = pen_tau;

