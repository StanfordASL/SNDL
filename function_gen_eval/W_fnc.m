function [W_p,W] = W_fnc(w_p, w_np, n, m)

%np' = (n-m)*((n-m)+1)/2
%n' = n*(n+1)/2

%w_p: Nc x np' matrix of Wp at Nc points encoded as vectors
%w_np: Nc x n' matrix of Wnp at Nc points encoded as vectors

%% Convert to matrices first

W_p = vec_to_sym(w_p,n-m);
W_np = vec_to_sym(w_np,n);

%% Assemble

W = W_np; %W_np is zero in top-left (n-m) block
W(1:(n-m), 1:(n-m), :) = W_p;
 


end