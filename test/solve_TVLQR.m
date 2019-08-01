function L_lqr = solve_TVLQR(t_ref,x_opt,u_opt,f_h,B_h,df_h,dB_h,n,m)

%x_dot = f(x) + B(x)*u
%delta x_dot = (dfdx + u*dB_dx)*delta_x + B(x)delta_u

Q = 4*eye(n);
R = eye(m);

%% Proceed with LQR

A = zeros(length(t_ref),n*n);
B = zeros(length(t_ref),n*m);

for i = 1:length(t_ref)
    A(i,:) = reshape(df_h(x_opt(i,:)'),1,n*n);
    dBdx = dB_h(x_opt(i,:)');
    for j = 1:m
        A(i,:) = A(i,:) + u_opt(i,j)*reshape(dBdx(:,:,j),1,n*n);
    end
    B(i,:) = reshape(B_h(x_opt(i,:)'),1,n*m);
end

t_bwd = t_ref; %backward time

% Riccati solution (backward)
[~,K]=ode45(@(tau,K)LQRfun(tau,K,t_ref,A,B,Q,R,n,m,t_ref(end)),t_bwd,reshape(Q,n*n,1));

K = flipud(K);
L_lqr=zeros(m,n,length(t_ref));

for i=1:length(t_ref)
    K_mat=reshape(K(i,:)',n,n);
    A_mat = reshape(A(i,:),n,n);
    B_mat = reshape(B(i,:),n,m);
    L_lqr(:,:,i) =-R\B_mat'*K_mat;
end

end

function dK_dt = LQRfun(tau,K,t_ref,A_ref,B_ref,Q,R,n,m,T)

A = reshape(interp1(t_ref,A_ref,T-tau),n,n);
B = reshape(interp1(t_ref,B_ref,T-tau),n,m);

K_mat = reshape(K,n,n);

dK_dt = reshape((Q-K_mat*B*(R\B')*K_mat +... 
    +K_mat*A+A'*K_mat),n*n,1);
end