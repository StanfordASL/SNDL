function [MP_state,MP_ctrl,converged] = ...
    compute_MP(Prob,act_p,n,m,N,L_e,XU_guess,t_ref,t_grid,x_con)

%% Solution guess
x_guess = interp1(t_ref,XU_guess(:,1:n),t_grid);
u_guess = interp1(t_ref,XU_guess(:,n+1:n+m),t_grid);

for k = 1:length(t_grid)
    for i = 1:n
        if x_guess(k,i) >= x_con(i,2)
            x_guess(k,i) = x_con(i,2)-(x_con(i,2)-x_con(i,1))*0.01;
        elseif x_guess(k,i) <= x_con(i,1)
            x_guess(k,i) = x_con(i,1)+(x_con(i,2)-x_con(i,1))*0.01;
        end
    end
end

x0 = reshape(x_guess',(N+1)*n,1);
u0 = reshape(u_guess',(N+1)*m,1);

Prob = modify_x_0(Prob,[x0;u0]);

%% Update constraint information
Prob.user.x_act = act_p;

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

end