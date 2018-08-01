function [W, dW] =  construct_metric(om_wp, theta_p, om_w_np, theta_np, O, n,m)


Phi_p_T =  @(X) (1/sqrt(O))*(kron(cos(X(:,1:n-m)*om_wp'),[1,0])+kron(sin(X(:,1:n-m)*om_wp'),[0,1]));
Phi_np_T = @(X) (1/sqrt(O))*(kron(cos(X*om_w_np'),[1,0])+kron(sin(X*om_w_np'),[0,1]));

W = @W_fun;

    function W = W_fun(X)
        W_p = vec_to_sym(Phi_p_T(X)*theta_p,n-m);
        W_np = vec_to_sym(Phi_np_T(X)*theta_np,n);
        
        W = W_np; %W_np is zero in top-left (n-m) block
        W(1:(n-m), 1:(n-m), :) = W_p;
    end

%%
dW = @Jacobian_W;

%dphi/dx_j (X):
Jw_p_j = @(j,X,N) (1/sqrt(O))* (kron(-sin(X(:,1:n-m)*om_wp').*kron(ones(N,1),om_wp(:,j)'),[1,0])+...
                                kron( cos(X(:,1:n-m)*om_wp').*kron(ones(N,1),om_wp(:,j)'),[0,1]));

Jw_np_j = @(j,X,N) (1/sqrt(O))* (kron(-sin(X*om_w_np').*kron(ones(N,1),om_w_np(:,j)'),[1,0])+...
                                 kron( cos(X*om_w_np').*kron(ones(N,1),om_w_np(:,j)'),[0,1]));


    function dWdx = Jacobian_W(X)
        
        N = size(X,1);
        
        dWdx = zeros(n,n,N,n);
        
        for j = 1:n-m
            %first construct dW_np part
            dWdx(:,:,:,j) = vec_to_sym(Jw_np_j(j,X,N)*theta_np,n); %zeros in top (n-m)x(n-m) block
            
            %add in dW_p part:
            dWdx(1:n-m,1:n-m,:,j) = vec_to_sym(Jw_p_j(j,X,N)*theta_p,n-m);
            
        end
        for j = n-m+1:n
            dWdx(:,:,:,j) = vec_to_sym(Jw_np_j(j,X,N)*theta_np,n);
        end
    end


end