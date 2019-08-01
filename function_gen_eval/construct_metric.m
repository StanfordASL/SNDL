function [W, dW_p] =  construct_metric(om_wp, theta_p, om_w_np, theta_np, O, n,m,w_const)


Phi_p_T =  @(X) (1/sqrt(O))*(kron(cos(X(:,1:n-m)*om_wp'),[1,0])+kron(sin(X(:,1:n-m)*om_wp'),[0,1]));
Phi_np_T = @(X) (1/sqrt(O))*(kron(cos(X*om_w_np'),[1,0])+kron(sin(X*om_w_np'),[0,1]));

W = @W_fun;

    function W = W_fun(X)
        
        N = size(X,1);
        
        W = vec_to_sym([Phi_p_T(X)*theta_p,Phi_np_T(X)*theta_np],n) +...
            repmat(w_const*eye(n),1,1,N);
        
%         W_p = vec_to_sym(Phi_p_T(X)*theta_p,n-m);
%         W_np = vec_to_sym(Phi_np_T(X)*theta_np,n);
%         W = W_np; %W_np is zero in top-left (n-m) block
%         W(1:(n-m), 1:(n-m), :) = W_p;
    end

%%
dW_p = @Jacobian_W;

%dphi/dx_j (X):
Jw_p_j = @(j,X,N) (1/sqrt(O))* (kron(-sin(X(:,1:n-m)*om_wp').*kron(ones(N,1),om_wp(:,j)'),[1,0])+...
                                kron( cos(X(:,1:n-m)*om_wp').*kron(ones(N,1),om_wp(:,j)'),[0,1]));


    function dWp_dx = Jacobian_W(X)
        
        N = size(X,1);
        
        dWp_dx = zeros(n-m,n-m,N,n-m);
        
        for j = 1:n-m
            %add in dW_p part:
            dWp_dx(:,:,:,j) = vec_to_sym(Jw_p_j(j,X,N)*theta_p,n-m);
        end
    end


end