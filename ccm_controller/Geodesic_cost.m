function J = Geodesic_cost(vars,w,n,K,W_fun,Phi,Phi_dot)

J = 0;

global GEO_X; %K+1 x n array of points
global GEO_MXDOT;

for k = 1:K+1 %0 ---> K
    GEO_X(k,:) = vars'*Phi(:,:,k)';
end

%Eval W at the GEO_X points
W = W_fun(GEO_X)+repmat(0.1*eye(n),1,1,K+1);

for k = 1:K+1 %0 ---> K
    x_dot_k = Phi_dot(:,:,k)*vars;
    
    M = W(:,:,k)\eye(n);
    
    M_xdot = M*x_dot_k;
    J = J + (1/2)*w(k)*(x_dot_k'*M_xdot); 
    
    GEO_MXDOT(:,k) = M_xdot;
end

return