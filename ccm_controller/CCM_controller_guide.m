%% Sample use of CCM controller

%% Setup offline

%parameters
%n: state-space dimension
n = 6;
%choose order of Chebyshev expansion for geodesics (each dimension)
geodesic_N = 2;

%contraction rate
lambda = 0.1;

%Functions
%W_fun: handle to vectorized function to compute W(x)
%dW_fun: handle to vectorized function to compute dW_xj
%f_fun: handle to compute dynamics function f(x)
%B_fun: handle to compute dynamics function B(x)

%initializes geodesic_MPC struct
[geo_Prob,geo_Ke,T_e,T_dot_e,geo_Aeq] = setup_geodesic_calc(n,geodesic_N,W_fun,dW_fun);

%choose solver
geo_solver = 'npsol';   

%initialize solution structure (used later for warm-start)
geo_warm = struct('sol',0,'result',[]);    

%% First time run (to generate solution structure)

% start_p = n x 1 start point of geodesic (nominal state)
% end_p = n x 1 end point of geodesic (actual state)

tic
[~, ~,J_opt,converged_geo,geo_result,geo_Prob] = ...
    compute_geodesic(geo_Prob,n,geodesic_N,...
            start_p,end_p,...
            T_e,T_dot_e,geo_Aeq,geo_warm,geo_solver);
toc;
disp('Geo dist: ');disp(converged_geo);
disp(sqrt(J_opt));
geo_Prob.CHECK = 1;
geo_warm.sol = 1;
geo_warm.result = geo_result;

%% Now ready for recursive online calls

%%%%% Repeat online %%%%%%

%@time t:
%geodesic from x_nom(t) -> x(t)
[X_e, X_dot_e,J_opt,converged_geo,geo_result,geo_Prob] = compute_geodesic(geo_Prob,...
            n,geodesic_N,x_nom,x_act,T_e,T_dot_e,geo_Aeq,geo_warm,geo_solver);
        
geo_warm.result = geo_result;

%compute control given u_nom
u_aux = compute_opt_aux(geo_Ke,X_e,X_dot_e,J_opt,W_fun,f,B,u_nom,lambda);

%implement:
u = u_nom + u_aux;

%%%%% Repeat online %%%%%%