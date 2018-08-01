%% Features

O_dyn = 8*n; %# random gaussian directions for dynamics
O_w = 6*n; %# random gaussian directions for W

D_dyn = 2*O_dyn*n;
D_w = 2*O_w;

%exclude these components of x in f(x)
n_dyn = 1:2;

%exclude these components of x in B(x)
n_dyn_B = 1:6; %B should be constant for PVTOL

%exclude these components of x in W(x)
n_W = 1:2;

%% Regularization constants

mu_f = 0.001;
mu_b = 0.001;
mu_s = 0.01;
mu_w = 0.0001;

%% Tolerance params

% slack tolerance for re-sample
tol = 0.05;

% extra tolerance for LMI constraints
eps_l = 0.1; %lambda
eps_wl = 0.1; %wl

% lower bounds on LMI constraints
delta_l = 0.1; %lambda
delta_wl = 0.1; %wl

%% Kernel definitions

kernel_f = struct('sigma',6,'L',eye(n));
kernel_b = struct('sigma',1,'L',blkdiag(zeros(n-m),eye(m)));
kernel_wp = struct('sigma',15.0);
kernel_wnp = struct('sigma',15.0);
