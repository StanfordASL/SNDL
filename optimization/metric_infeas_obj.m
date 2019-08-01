function [obj,prob] = metric_infeas_obj(vars,prob,theta_p_prev,theta_np_prev,Zw_p,Zw_np,Jw_p,f,dfdX_p,Xc_i,D,n,m,Nc,...
                            n_p,w_lower,delta_wl,eps_wl,lambda,eps_l,pen_t)

%min violation to get strictly feasible point

%%


% phase 1 cost obj regularization
mu_w_phase1 = prob.mu_w_phase1;


%% Extract variables

theta_p = vars(1:D*n_p);
theta_np = vars(D*n_p+1:end);

Theta_p = reshape(theta_p,D,[]);
Theta_np = reshape(theta_np,D,[]);

Theta_p_prev = reshape(theta_p_prev,D,[]);
Theta_np_prev = reshape(theta_np_prev,D,[]);

%% Compute W_p, W, partial_Wp

[W_p,W] = W_fnc(Zw_p*Theta_p,Zw_np*Theta_np,n,m,w_lower);
partial_Wp = dWp_fnc(Jw_p,Theta_p,Nc,n,m);

dWp_f_all = zeros((n-m),(n-m),Nc);
for j = 1:n-m
    dWp_f_all = dWp_f_all + partial_Wp(:,:,:,j).*reshape(f((Xc_i-1)*n+j),1,1,Nc);
end

%% Constraints

MAX_F = zeros(Nc,1);
PENS_F = zeros(Nc,1);
EIG_V = zeros(n-m,n-m,Nc);
EIG_L = zeros(Nc,n-m);

MAX_W = zeros(Nc,1);
PENS_W = zeros(Nc,1);
EIG_V_W = zeros(n,n,Nc);
EIG_L_W = zeros(Nc,n);

z_rand_W = mvnrnd(zeros(n,1),eye(n),Nc);
z_rand_F = mvnrnd(zeros(n-m,1),eye(n-m),Nc);


for k = 1:Nc
    
    %% eigenvalues of W
    
    [V_W, E_W] = eig((delta_wl+eps_wl)*eye(n) - W(:,:,k));
    [Es_W,order_w] = sort(diag(E_W));
    
    if (sum(Es_W(1:n-1)==Es_W(n))>0)   
        [V_W, E_W] = eig((delta_wl+eps_wl)*eye(n) - W(:,:,k) + (0.00001/n)*(z_rand_W(k,:)'*z_rand_W(k,:)));
        [Es_W,order_w] = sort(diag(E_W));
    end
    
    %arrange in ascending order
    Vs_W = V_W(:,order_w);
    
    %compute violation
    MAX_W(k) = Es_W(n);
    
    %compute penalty
    PENS_W(k) = (MAX_W(k) > pen_t) * (MAX_W(k) - pen_t)^3;
        
	%save
    EIG_L_W(k,:) = Es_W';
    EIG_V_W(:,:,k) = Vs_W;
    
    %% eigenvalues of F
    
    dWp_f = dWp_f_all(:,:,k);
    
    F = -dWp_f + dfdX_p(:,:,k)*W(:,1:n-m,k) + W(1:n-m,:,k)*dfdX_p(:,:,k)' + 2*(lambda+eps_l)*W_p(:,:,k);
    
    %eigvalues
    [V,E] = eig(F);
    [Es,order] = sort(diag(E));
    
    if (sum(Es(1:n-m-1)==Es(n-m))>0)  
        [V,E] = eig(F + (0.00001/(n-m))*(z_rand_F(k,:)'*z_rand_F(k,:)));
        [Es,order] = sort(diag(E));
    end
	
	%sort in ascending order
    Vs = V(:,order);
    
    %largest eig value
    MAX_F(k) = Es(n-m);
    
    %penalty
    PENS_F(k) =  max(MAX_F(k)-pen_t,0.0); 
    
    EIG_L(k,:) = Es';
    
    EIG_V(:,:,k) = Vs;
          
end

%% Combine

obj = sum(PENS_W) + mu_w_phase1*sum(PENS_F) + ...
                    mu_w_phase1*trace((Theta_p-Theta_p_prev)'*(Theta_p-Theta_p_prev)) +...
                    mu_w_phase1*trace((Theta_np-Theta_np_prev)'*(Theta_np-Theta_np_prev));

%% Prob struct output

prob.MAX_F = MAX_F;
prob.PENS_F = PENS_F;
prob.EIG_L = EIG_L;
prob.EIG_V = EIG_V;

prob.MAX_W = MAX_W;
prob.PENS_W = PENS_W;
prob.EIG_V_W = EIG_V_W;
prob.EIG_L_W = EIG_L_W;


end