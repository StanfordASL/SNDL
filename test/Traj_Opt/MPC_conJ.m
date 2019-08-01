function conJ = MPC_conJ(xu,Prob,n,m,N,P,D,df,dB,Tp)
%Dynamics, and terminal

global B_full

conJ = zeros(n*(N+1)+2,(n+m)*(N+1));

%% Dynamics constraints

conJ(1:n*(N+1),1:n*(N+1)) = (2/Tp)*D - ...
            df_all(df,xu(1:n*(N+1)),n,N) - ...
            dB_all(dB,xu(1:n*(N+1)),xu(n*(N+1)+1:end),n,m,N);

conJ(1:n*(N+1),n*(N+1)+1:end) = -B_full;

%% Initial state constraint
x_init = Prob.user.x_act;
conJ(n*(N+1)+1,1:n) = 2*(xu(1:n)-x_init)';

%% Terminal constraint

x_eq = Prob.user.x_eq;
conJ(n*(N+1)+2,n*N+1:n*(N+1)) = 2*(xu(n*N+1:n*(N+1))-x_eq)'*P;

end

