function [W_p,W] = W_fnc(w_p, w_np, n, m, w_const)

%np' = (n-m)*((n-m)+1)/2
%n' = n*(n+1)/2

%w_p: Nc x np' matrix of Wp at Nc points encoded as vectors
%w_np: Nc x n' matrix of Wnp at Nc points encoded as vectors

%% Convert to matrices first

N = size(w_p,1);

W = vec_to_sym([w_p,w_np],n) + repmat(w_const*eye(n),1,1,N);

W_p = W(1:(n-m),1:(n-m),:);
 


end