function c = MP_con(xu,Prob,n,m,N,P,D,f,B,Tp,x_eq)
%Dynamics, and terminal

%actual initial state
x_init = Prob.user.x_act;

c = zeros(n*(N+1)+2,1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Tp)*D*xu(1:n*(N+1)) -...
                dyn_all(f,B,xu(1:n*(N+1)),xu(n*(N+1)+1:end),n,m,N);

%% Initial state constraint

c(n*(N+1)+1) = (xu(1:n)-x_init)'*(xu(1:n)-x_init);

%% Terminal constraint
c(n*(N+1)+2) = (xu(n*N+1:n*(N+1))-x_eq)'*P*(xu(n*N+1:n*(N+1))-x_eq);

end

