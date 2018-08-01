function GradObj = Geodesic_grad(~,w,K,N,n,Ti,dW_fun,Phi_dot)

global GEO_X; %K+1 x n
global GEO_MXDOT; %n x K+1

GradObj = zeros(n*(N+1),1);

%eval dW at GEO_X points
%dW(:,:,k,j) = dW_xj (x_k)
dW = dW_fun(GEO_X); 

for k = 1:K+1 %0 ---> K
    M_xdot = GEO_MXDOT(:,k);
    
    GradObj = GradObj + w(k)*Phi_dot(:,:,k)'*M_xdot;
    
%     W_dx = dW(GEO_X(:,k));
    
    for j = 1:n
        GradObj = GradObj+...
           -(1/2)*w(k)*(M_xdot'*dW(:,:,k,j)*M_xdot)*Ti(:,k,j);
    end
end


return