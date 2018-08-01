function [J_opt,converged,Prob,controller,u_aux,ccm_ok,lambda_act,d_est] = compute_CCM_controller(geo_Prob,info,x_nom,u_nom,x,...
                                                                            ok_prev,lambda_prev,d_prev)


%compute CCM feedback controller
%outputs:
%J_opt: geodesic energy E
%converged: geodesic convergence flag
%controller: updated controller struct
%u_aux: feedback
%ccm_ok: all computations ok
%lambda_act: actual lambda used for computing u_ccm
%d_est: estimated mismatch between true model and CCM learned model

ccm_thresh = info.E_thresh;

%% Compute geodesic

[X_e,X_dot_e,J_opt,converged,geo_result,Prob] = ...
    compute_geodesic(geo_Prob,info.n,info.N,x_nom,x,info.T_e,info.T_dot_e,info.Aeq,info.warm,info.solver);

%% Update struct for next warm start

% fprintf('ok: %d, E: %.4f\n',converged, J_opt);

controller = info;
if (converged == 0 || converged == 1) 
    if (ok_prev) %previous comp was valid
        %get previous E
        E_prev = info.warm.E;
        if E_prev < ccm_thresh %i.e. last comp valid but did not use CCM
            lambda_prev = 0;
        end            
    else
        E_prev = J_opt;
        lambda_prev = 0;
    end
    controller.warm.sol = 1;
    controller.warm.result = geo_result;
    controller.warm.E = J_opt;
else
    controller.warm.sol = 0;
    u_aux = zeros(info.m,1);
    lambda_act = -1;
    ccm_ok = 0;
    d_est = 0;
    return;
end

%% Decide whether to do CCM control

if (J_opt < ccm_thresh)
    u_aux = zeros(info.m,1);
    ccm_ok = 1; %to allow recovery of E_prev
    lambda_act = 0; %no ccm control
    d_est = 0;
    return;
end

%% Adaptively choose lambda based on disturbance estimate

[lambda_act, d_est] = check_lambda(info.n,info.m,X_e,X_dot_e,info.f,info.B,info.df,info.W,info.dW,J_opt,E_prev,d_prev,lambda_prev,info.lambda,info.dt,ccm_thresh);

if lambda_act < 0 %ccm controller not valid
    u_aux = zeros(info.m,1);
    ccm_ok = 0;
    d_est = 0;
    return;
else
    ccm_ok = 1;
end

%% Compute control

%only if successful geodesic solve and lambda_act >= 0
u_aux = compute_opt_aux(info.Ke,X_e,X_dot_e,J_opt,info.W,info.f,info.B,u_nom,lambda_act);

end