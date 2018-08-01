function dW_p = dWp_fnc(Jw_p,theta_p,Nc,n,m)

%Jw_p: Nc x D x n-m array: Jw_p(:,:,j) = partial_j phi_p^T (X)
%dW_p: n-m x n-m x Nc x n-m array: dW(:,:,:,j) = partial_j W_p(X)

dW_p = repmat(theta_p(1,1),n-m,n-m,Nc,n-m);

for j = 1:n-m
    dW_p(:,:,:,j) = vec_to_sym(Jw_p(:,:,j)*theta_p,n-m);
end



end