function [Phi_T,W] = Phi_b_T(X,O,kernel,n_B)

%%%%%% Info %%%%%%
%phi(x): D x n random feature matrix for:
%approximation: b_j(x) = phi_T(x) * beta_j

%%%%%% Inputs %%%%%%
%X: N x n matrix of state points
%O: #gaussian directions
%kernel: kernel params
%n_B: ignore these components for B

%%%%%% Outputs %%%%%%
%Phi_T: Nn x D matrix: b_j(X) = Phi_T(X) * beta_j
%W: sampled frequencies for random features

%% constants
n = size(X,2);

%% Feature matrix
% Defined by: Stack: (cos(om'*x); sin(om'*x))*L(om)
% om ~ sampled from rho

%% Gaussian separable
%positive matrix part of kernel
L = kernel.L;

%gaussian kernel variance
sigma = kernel.sigma;

% rho ~ N(0, (1/sigma^2)*eye(n))
% sample O points
W = mvnrnd(zeros(n,1),(1/sigma^2)*eye(n),O); %O x n
W(:,n_B) = 0;

XW_ = X*W'; % N x O

Phi_T = (1/sqrt(O)) * ( kron ( kron(cos(XW_),[1,0]) + kron(sin(XW_),[0,1]) , L' ) );

end