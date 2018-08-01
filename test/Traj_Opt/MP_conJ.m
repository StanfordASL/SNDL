function conJ = MP_conJ(xu,Prob)%,D,n,m,N,df,B_full,Tp,P)
%Dynamics, and terminal
n = Prob.user.n;
N = Prob.user.N;
m = Prob.user.m;

conJ = zeros(n*(N+1)+2,(n+m)*(N+1));

%% Dynamics constraints

conJ(1:n*(N+1),1:n*(N+1)) = (2/Prob.user.Tp)*Prob.user.D - ...
            df_all(Prob.user.df,xu(1:n*(N+1)),n,N);

conJ(1:n*(N+1),n*(N+1)+1:end) = -Prob.user.B_full;

%% Initial state constraint

conJ(n*(N+1)+1,1:n) = 2*(xu(1:n)-Prob.user.x_act)';

%% Terminal constraint
conJ(n*(N+1)+2,n*N+1:n*(N+1)) = 2*((xu(n*N+1:n*(N+1))-Prob.user.x_eq))';

end

