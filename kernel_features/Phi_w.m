function [Phi_p_T_all,Phi_p_T,Phi_np_T_all,Phi_np_T,Jw_p_all,Jw_p,Wp,Wnp] = Phi_w(n_W,X,Xc_i,...
                                        O,kernel_p,kernel_np,m,flag_sample_om,Wp,Wnp)

%%%%%% Info %%%%%%
%phi_p(x): D x n random feature matrix for:
%approximation: Wp_ij(x) = phi_p(x) * theta_ij

%phi_np(x): D x n random feature matrix for:
%approximation: Wnp_ij(x) = phi_np(x) * theta_ij

%W(x) = [Wp(x)] + Wnp(x)
%Wp : (n-m) x (n-m), [ ]->project up to n x n
%Wnp: n x n, with (n-m) x (n-m) top-left block all zeros

%%%%%% Inputs %%%%%%
%n_W: ignore these components for W
%X: N x n matrix of state points
%O: #gaussian directions
%Xc_i: constraint set indices in X
%kernel_p, kernel_np: kernel params
%flag_sample_om: prevent re-sampling of gaussian directions
%Wp, Wnp: existing gaussian samples (if they exist)
%m: #inputs

%%%%%% Outputs %%%%%%
%Wp, Wnp: gaussian samples
%Phi_p_T_all: N x D matrix:  Phi_p_T(X)
%Phi_np_T: N x D matrix: Phi_np_T(X)
%Phi_P_T, Phi_np_T: Nc x D subsample for current Xc set
%Jw_p_all: N x D x n-m array: Jw_p(:,:,j) = D_j phi_p^T (X) 
%Jw_p: Nc x D x n-m subsample for current Xc set

%% constants
n = size(X,2);
N = size(X,1);

%% Feature matrix
% Defined by: Stack: (cos(om'*x); sin(om'*x))
% om ~ sampled from rho

%% Gaussian 
%gaussian kernel variance
sigma_p = kernel_p.sigma;
sigma_np = kernel_np.sigma;

% rho ~ N(0, (1/sigma^2)*eye)
% sample O points
if (flag_sample_om)
    Wp = mvnrnd(zeros(n-m,1),(1/sigma_p^2)*eye(n-m),O); %O x n-m
    Wnp  = mvnrnd(zeros(n,1), (1/sigma_np^2)*eye(n),O); %O x n
    
    Wp(:,n_W) = 0;
    Wnp(:,n_W) = 0;
end
D = 2*O;

XW_p = X(:,1:n-m)*Wp'; % N x O
XW_np = X*Wnp'; %N x O

Phi_p_T_all =  (1/sqrt(O))*( kron(cos(XW_p),[1,0]) + kron(sin(XW_p),[0,1]) );
Phi_np_T_all = (1/sqrt(O))*( kron(cos(XW_np),[1,0]) + kron(sin(XW_np),[0,1]) );

Phi_p_T = Phi_p_T_all(Xc_i,:);
Phi_np_T = Phi_np_T_all(Xc_i,:);

%% Gradient matrix

Jw_p_all = zeros(N,D,n-m);

for j = 1:n-m
    
    Jw_p_all(:,:,j) = (1/sqrt(O))* (kron(-sin(XW_p).*kron(ones(N,1),Wp(:,j)'),[1,0])+...
                                    kron( cos(XW_p).*kron(ones(N,1),Wp(:,j)'),[0,1]));
end

Jw_p = Jw_p_all(Xc_i,:,:);

end