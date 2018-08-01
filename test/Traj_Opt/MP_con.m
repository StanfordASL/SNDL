function c = MP_con(xu,Prob)
%Dynamics, and terminal

n = Prob.user.n;
N = Prob.user.N;

c = zeros(n*(N+1)+2,1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Prob.user.Tp)*Prob.user.D*xu(1:n*(N+1)) -...
    (f_all(Prob.user.f,xu(1:n*(N+1)),n,N) + Prob.user.B_full*xu(n*(N+1)+1:end));

%% Initial state constraint

c(n*(N+1)+1) = (xu(1:n)-Prob.user.x_act)'*eye(n)*(xu(1:n)-Prob.user.x_act);

%% Terminal constraint
c(n*(N+1)+2) = (xu(n*N+1:n*(N+1))-Prob.user.x_eq)'*(xu(n*N+1:n*(N+1))-Prob.user.x_eq);

end

