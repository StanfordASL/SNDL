function u_aux = compute_opt_aux(K_e,X,X_dot,E,...
                            W_fun,f_fun,B_fun,u_nom,lambda)

%K_e: number of points in geodesic trace X
%X: n x K_e matrix: [x_1....x_{Ke}] geodesic points
%X_dot: n x K_e matrix: [x_dot_1...x_dot_{Ke}] geodesic curve velocities

%%

W = W_fun(X');

k = K_e;
A = 2*X_dot(:,k)'*(W(:,:,k)\B_fun(X(:,k)));

% l_b = -Inf;
u_b = -2*lambda*E + 2*X_dot(:,1)'*(W(:,:,1)\(f_fun(X(:,1)) + B_fun(X(:,1))*u_nom)) -...
                    2*X_dot(:,k)'*(W(:,:,k)\(f_fun(X(:,k)) + B_fun(X(:,k))*u_nom));

a = -u_b;
b = A';
if (norm(b)==0) || (a <= 0) 
    u_aux = zeros(length(b),1);
else
    u_aux = -(a/(b'*b))*b;
end

end



