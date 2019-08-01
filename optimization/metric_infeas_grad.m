function [J,H] = metric_infeas_grad(vars,prob,theta_p_prev,theta_np_prev,Nc,dW_dtheta,dF_dtheta,n,m,n_vars,pen_t)

%min violation to get strictly feasible point


% phase 1 cost obj regularization
mu_w_phase1 = prob.mu_w_phase1;

%%

MAX_F = prob.MAX_F;
PENS_F = prob.PENS_F;
EIG_V = prob.EIG_V;
EIG_L = prob.EIG_L;

MAX_W = prob.MAX_W;
PENS_W = prob.PENS_W;
EIG_V_W = prob.EIG_V_W;
EIG_L_W = prob.EIG_L_W;

I_W = eye(n);
I_F = eye(n-m);

%% Jacobian & Hessian from regularization

J = 2 * mu_w_phase1 * (vars-[theta_p_prev;theta_np_prev]);
H = 2 * mu_w_phase1 * eye(n_vars);

%% Gradients from Contraints

J_constraints = zeros(n_vars,Nc);
H_constraints = zeros(n_vars,n_vars,Nc);

for k = 1:Nc
    
	%% Jacobian & Hessian from W Constraints
    if PENS_W(k) > 0
        
        %Jacobian
        
        v = EIG_V_W(:,n,k);

        dW_dtheta_v = kron(v',I_W)* (-dW_dtheta(:,:,k));
        d_eig_W_dtheta = v'*(dW_dtheta_v);

        pen_j_W = 3.0 * (MAX_W(k) - pen_t)^2;

        J_constraints(:,k) = pen_j_W * d_eig_W_dtheta';

        %Hessian

        l = EIG_L_W(k,n);

        Quad_Prod_W =  dW_dtheta_v'*EIG_V_W(:,1:n-1,k);
        H_W = 2*pen_j_W* (Quad_Prod_W *diag( 1./(l-EIG_L_W(k,1:n-1)) )*Quad_Prod_W');

        pen_jj_W = 6.0 * (MAX_W(k) - pen_t);
        H_W = H_W + pen_jj_W * (d_eig_W_dtheta' * d_eig_W_dtheta);

        H_constraints(:,:,k) = H_W;
        
    end
	
	%% Jacobian & Hessian from F Constraints
	
	if PENS_F(k)>0
		
		%Jacobian
		
		v = EIG_V(:,n-m,k);
        
        dF_dtheta_v = kron(v',I_F)*dF_dtheta(:,:,k);
        d_eig_F_dtheta = v'*dF_dtheta_v;        

        
        pen_j =  1.0;
        
        J_constraints(:,k) = J_constraints(:,k) +  mu_w_phase1 * pen_j * (d_eig_F_dtheta');
		
		
		%Hessian
		l = EIG_L(k,n-m);
        
        Quad_Prod_F =  dF_dtheta_v'*EIG_V(:,1:n-m-1,k);
        H_F = 2*pen_j* (Quad_Prod_F *diag( 1./(l-EIG_L(k,1:n-m-1)) )*Quad_Prod_F');
        
        
        H_constraints(:,:,k) = H_constraints(:,:,k) +  mu_w_phase1 * H_F;
	
	end
	
	
    
end


%% assemble

J = J + sum(J_constraints,2);

H = H + sum(H_constraints,3);

end