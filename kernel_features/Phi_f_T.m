function [Phi_T,Jf_all,Jf,W] = Phi_f_T(X,O,Xc_i,kernel,m,n_dyn)

%%%%%% Info %%%%%%
%phi(x): D x n random feature matrix for:
%approximation: f(x) = phi_T(x) * alpha

%%%%%% Inputs %%%%%%
%X: N x n matrix of state points
%O: #gaussian directions
%Xc_i: constraint set indices in X
%kernel: kernel params
%m: #inputs
%n_dyn: ignore these components for f

%%%%%% Outputs %%%%%%
%Phi_T: Nn x D matrix: f(X) = Phi_T(X) * alpha
%Jf_all: D x n x N x (n-m) array: Jf_all(:;,:,i,j) = D phi_j(x_i)
%Jf: D x n x Nc x (n-m) array: sub-sample of Jf_all for current Xc
%W: sampled frequencies for random features

%% constants
n = size(X,2);
N = size(X,1);

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
W(:,n_dyn) = 0;
D = 2*O*n;

XW_ = X*W'; % N x O

Phi_T = (1/sqrt(O)) * ( kron ( kron(cos(XW_),[1,0]) + kron(sin(XW_),[0,1]) , L' ) );

%% Gradient matrix

Jphi_j = @(x) kron(-diag(sin(W*x))*W,[1;0]) + kron(diag(cos(W*x))*W,[0;1]);

Jf_all = zeros(D,n,N,n-m);

for i = 1:N
    x = X(i,:)';
    for j = 1:n-m
        Jf_all(:,:,i,j) = (1/sqrt(O))*kron(Jphi_j(x),L(:,j));
    end
end

Jf = Jf_all(:,:,Xc_i,:);


end