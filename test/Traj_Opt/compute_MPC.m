function [MP_state,MP_ctrl,converged,mpc_result] = ...
    compute_MPC(Prob,act_p,n,m,N,L_e,XU_nom,t_nom,t_grid,x_con,u_con,mpc_result)


if (~mpc_result.sol)
    %% Solution guess
    x_guess = interp1(t_nom,XU_nom(:,1:n),t_grid);
    u_guess = interp1(t_nom,XU_nom(:,n+1:n+m),t_grid);
    
    for k = 1:length(t_grid)
        for i = 1:n
            if x_guess(k,i) >= x_con(i,2)
                x_guess(k,i) = x_con(i,2)-(x_con(i,2)-x_con(i,1))*0.05;
            elseif x_guess(k,i) <= x_con(i,1)
                x_guess(k,i) = x_con(i,1)+(x_con(i,2)-x_con(i,1))*0.05;
            end
        end
        for j = 1:m
            if u_guess(k,j) >= u_con(j,2)
                u_guess(k,j) = u_con(j,2)-(u_con(j,2)-u_con(j,1))*0.05;
            elseif u_guess(k,j) <= u_con(j,1)
                u_guess(k,j) = u_con(j,1)+(u_con(j,2)-u_con(j,1))*0.05;
            end
        end
    end

    x0 = reshape(x_guess',(N+1)*n,1);
    u0 = reshape(u_guess',(N+1)*m,1);

    Prob = modify_x_0(Prob,[x0;u0]);

else
    Prob = WarmDefSOL('snopt',Prob,mpc_result.result);
end

%% Update constraint information
x_nom = interp1(t_nom,XU_nom(:,1:n),t_grid);
u_nom = interp1(t_nom,XU_nom(:,n+1:n+m),t_grid);
x_nom_all = reshape(x_nom',(N+1)*n,1);
u_nom_all = reshape(u_nom',(N+1)*m,1);
xu_nom = [x_nom_all;u_nom_all];

Prob.user.x_act = act_p;
Prob.user.xu_nom = xu_nom;
Prob.user.x_eq = x_nom(end,:)';

Prob = ProbCheck(Prob,'snopt');

%% Solve

Result = snoptTL(Prob);
converged = Result.Inform; %GOOD: {1,2,3}

%% Compute trajectories

MP_state = zeros(size(L_e,2),n);
for i = 1:n
    c = Result.x_k(i:n:n*(N+1)-(n-i))';
    MP_state(:,i) = (c*L_e)';
end


MP_ctrl = zeros(size(L_e,2),m);
for j = 1:m
    c = Result.x_k(n*(N+1)+j:m:end-(m-j))';
    MP_ctrl(:,j) = (c*L_e)';
end

MP_ctrl = MP_ctrl - XU_nom(:,n+1:n+m);

%% If solve successfully, update result
if (converged==1)
    mpc_result.sol = 1;
    mpc_result.result = Result;
end

end