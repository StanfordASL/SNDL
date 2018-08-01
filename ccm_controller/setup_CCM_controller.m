function [geo_Prob,cntrl_info] = setup_CCM_controller(n,m,geodesic_N,...
                                                      x_lim_n,x_lim,...
                                                      W_fun,dW_fun,f_fun,B_fun,df_fun,...
                                                      lambda,x_0,u_0,dt,E_thresh)

% setup CCM controller 

%inputs:
%n,m: state and control dimension
%geodesic_N: order of Chebyshev approximation for geodesic
%x_lim_n,x_lim: indices and bounds for absolute limits on geodesic
%W_fun, dW_fun: function handles for W(x) and dW_dx
%f_fun,B_fun,df_fun: function handles for f(x), B(x), dfdx(x)
%lambda: contraciton rate
%x_0, u_0: test initial nominal state and control
%dt: zero-order-hold timestep for controller
%E_thresh: threshold for switching to CCM controller

%outputs: 
% geo_Prob: geodesic problem structure (Tomlab);
% cntrl_inf: controller struct

%% Setup geodesic solver

[geo_Prob,geo_Ke,T_e,T_dot_e,geo_Aeq] = setup_geodesic_calc(n,x_lim_n,x_lim,geodesic_N,W_fun,dW_fun);

%choose solver
geo_solver = 'npsol';   

%initialize solution structure (used later for warm-start)
geo_warm = struct('sol',0,'result',[],'E',0);

%define controller struct (maybe not most efficient data type in matlab)
cntrl_info.n = n;
cntrl_info.m = m;
cntrl_info.Ke = geo_Ke;
cntrl_info.T_e = T_e;
cntrl_info.T_dot_e = T_dot_e;
cntrl_info.Aeq = geo_Aeq;
cntrl_info.N = geodesic_N;
cntrl_info.solver = geo_solver;
cntrl_info.warm = geo_warm;
cntrl_info.lambda = lambda;
cntrl_info.W = W_fun;
cntrl_info.dW = dW_fun;
cntrl_info.f = f_fun;
cntrl_info.df = df_fun;
cntrl_info.B = B_fun;
cntrl_info.dt = dt;
cntrl_info.E_thresh = E_thresh;

%% First time run (to generate solution structure)

% start_p = n x 1 start point of geodesic (nominal state)
% end_p = n x 1 end point of geodesic (actual state)

x = x_0 + 0.1*randn(n,1);

tic
[J_opt,converged_geo,geo_Prob,cntrl_info,u_aux,ccm_ok,lamb,d]  = compute_CCM_controller(geo_Prob,cntrl_info,x_0,u_0,x,0,lambda,0);
toc;

fprintf('Geodesic:%d (good: 0, 1, 6), Geodesic dist: %.3f, u_fb: \n',converged_geo,sqrt(J_opt));
disp(u_aux);
geo_Prob.CHECK = 1;


end