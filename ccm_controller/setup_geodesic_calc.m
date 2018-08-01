function [geo_Prob,K_e,T_e,T_dot_e,Aeq] = ...
    setup_geodesic_calc(n,x_lim_n,x_lim,N,W,dW)
%n: state-space dimension
%N: Chebyshev polynomial order
%W, dW: W and dW matrix functions

%geo_Prob: Tomlab problem structure
%T_e: evaluation matrix for geodesic curve
%T_dot_e: evaluation matrix for geodesic curve velocity

K = 2*N; %number of computation points for cost,gradient,etc
K_e = 5; %number of evaluation points of final geodesic (beginning,middle,end)

%Optimization variables: chebyshev coefficients for each dimension for geodesic
%{c_0^1...c_N^1},...,{c_0^1...c_N^n}

%% Obtain Chebyschev Pseudospectral Numerics

%CGL points and quadrature weights
[t,w] = clencurt(K);

[t_e,~] = clencurt(K_e-1);

%compute Cheby basis at all points
[T, T_dot] = ...
    compute_cheby(K,N,t);

[T_e, T_dot_e] = ...
    compute_cheby(K_e-1,N,t_e);

% Chebyshev polynomial method
[phi_start, ~] = compute_cheby(0,N,-1);
[phi_end, ~] = compute_cheby(0,N,1);
A_start = kron(eye(n),phi_start');
A_end = kron(eye(n),phi_end');

%construct limit constraints
In = eye(n);
Lim_n = In(x_lim_n,:); %which components are limited
A_lim = [];
for k = 1:K+1
    A_lim = [A_lim;
             kron(Lim_n,T(:,k)')];
end
b_lim = repmat(x_lim,K+1,1);


%for initial and final points of geodesic
Aeq = sparse([A_start; %n
              A_end;   %n
              A_lim]); %length(b_lim)
          
b_L = [zeros(2*n,1);-b_lim];
b_U = [zeros(2*n,1); b_lim];

%use to evaluate x_k, x_dot_k
Phi = zeros(n,n*(N+1),K+1);
Phi_dot = zeros(n,n*(N+1),K+1);
for k = 1:K+1
    Phi(:,:,k) = kron(eye(n),T(:,k)');
    Phi_dot(:,:,k) = 2*kron(eye(n),T_dot(:,k)');
end

%use to evaluate cost derivative
Ti = zeros(n*(N+1),K+1,n);
In = eye(n);
for i= 1:n
    Ti(:,:,i) = kron(In(:,i),T);
end

%Yes global vars are bad, but useful here for memory sharing between
%Tomlab's subroutines
global  GEO_X; 
GEO_X = zeros(K+1,n);

global  GEO_MXDOT; 
GEO_MXDOT = zeros(n,K+1);

% Cost function
geo_cost_fnc =  @(vars) Geodesic_cost(vars,w,n,K,W,Phi,Phi_dot);

%Gradient function
geo_grad_fnc = @(vars) Geodesic_grad(vars,w,K,N,n,Ti,dW,Phi_dot);

Name = 'Geodesic Problem';
geo_Prob = conAssign(geo_cost_fnc,geo_grad_fnc,[],[],...
                  [],[],Name,zeros(n*(N+1),1),[],0,...
                  Aeq,b_L,b_U,[],[],[],[],[],[]);
              
% function Prob = conAssign(f, g, H, HessPattern, x_L, x_U, Name, x_0, ...
%                           pSepFunc, fLowBnd, ...
%                           A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ...
%                           x_min, x_max, f_opt, x_opt);

end