function [pos_coef,vel_coef,acc_coef,jer_coef,T_seg,coeff] = smooth_2D_path(path,n_traj,varargin)
%path: FMT* (or waypoint) path

%assume all same length paths
% n_paths = n_traj/n_demo;
N_nodes = size(path,1);
N_seg = N_nodes-1;

if nargin > 2
    v_init = varargin{1};
else
    v_init = zeros(n_traj,2);
end

%% Setup polynomial structure

%order
p_order = 10;

%coefficient return for pos,vel,acc,jerk
pos_coef = @(T) [1, T, T^2, T^3, T^4, T^5, T^6, T^7, T^8, T^9];
vel_coef = @(T) [0, 1, 2*T, 3*T^2, 4*T^3, 5*T^4, 6*T^5, 7*T^6, 8*T^7, 9*T^8];
acc_coef = @(T) [0, 0, 2, 6*T, 12*T^2, 20*T^3, 30*T^4, 42*T^5, 56*T^6, 72*T^7];
jer_coef = @(T) [0, 0, 0, 6, 24*T, 60*T^2, 120*T^3, 210*T^4,336*T^5,504*T^6];
snap_coef = @(T) [0,0,0,0,24,120*T,360*T^2,840*T^3,1680*T^4,3024*T^5];

%% Get Min Snap Q matrix

p = sym('p',[p_order,1]);
syms T
snap_poly_sq = (snap_coef(T)*p).^2;
J = int(snap_poly_sq,T);

Q = sym('Q',[p_order,p_order]);
Q(1:4,:) = zeros(4,p_order);
for i = 5:p_order
    for j = i:p_order
        Q(i,j) = diff(diff(J,p(i)),p(j))/2;
    end
end
for i = 5:p_order
    for j = 1:i-1
        Q(i,j) = Q(j,i);
    end
end
Q_fun = matlabFunction(Q); %symbolic function

%% Get Min Time Gram matrix

% Tv = sym('Tv',[N_seg,1]);
% Q_full = Q_fun(Tv(1));
% for ns = 2:N_seg
%     Q_full = blkdiag(Q_full,Q_fun(Tv(ns)));
% end
% A = get_constraints(N_seg,Tv,p_order,pos_coef,vel_coef,acc_coef,jer_coef,snap_coef);
% A_inv = inv(A);
% G = A_inv.'*Q_full*A_inv;
%
% %now convert to matlab functions
% Q_fun = matlabFunction(Q);
% G_T_full = cell(N_seg,1);
% G_T_fun = cell(N_seg,1);
% for ns = 1:N_seg
%     G_T_full{ns} = diff(G,Tv(ns));
%     G_T_fun{ns} = matlabFunction(G_T_full{ns});
% end

%% Minimize snap

coeff = cell(n_traj,1);
T_seg = cell(n_traj,1);
for i = 1:n_traj    
    %guess some segment times (will give slightly different snap trajectories centered on same path)
    T_seg{i} = zeros(N_seg,1);
    for ns = 1:N_seg
        T_seg{i}(ns) = (0.85+2.0*rand(1))*norm(path(ns+1,:)-path(ns,:))/2.0;
    end
    
    v0 = v_init(i,:);
    
    d_wp = zeros(2*N_seg,2);
    for ns = 1:N_seg
        d_wp(1+(ns-1)*2,:) = path(ns,:);
        d_wp(2*ns,:) = path(ns+1,:);
    end
    d_F = [d_wp;
           v0;
           zeros(2+2+4*(N_seg-1),2)];
    s_F = size(d_F,1);
    s_P = N_seg*p_order - s_F;
    
    %% Extract optimal solution
    
    [d_P, A_inv, ~] = get_dP(d_F,s_F,s_P,N_seg,T_seg{i},p_order,Q_fun,pos_coef,vel_coef,acc_coef,jer_coef,snap_coef);
    coeff{i} = A_inv*[d_F;d_P];
    
end

%% Compute polynomials

% poly_seg = cell(N_seg,1);
% dt = 0.01;
% for ns = 1:N_seg
%     len = floor(T_seg(ns)/dt)+1;
%     poly_seg{ns} = zeros(len,2);
%     for i = 1:len
%         t = (i-1)*dt;
%         poly_seg{ns}(i,:) = pos_coef(t)*coeff(1+(ns-1)*p_order:ns*p_order,:);
%     end
% end

%% Plot
% close all;
% figure(2)
% title('waypoint + polyspline path');
% hold on
% for i = 1:length(path)-1
%     line([path(j,1), path(j+1,1)],[path(j,2),path(j+1,2)],...
%             'Color', '-bx', 'LineWidth', 2);
% end
% for ns = 1:N_seg
%     plot(poly_seg{ns}(:,1),poly_seg{ns}(:,2),'r-','linewidth',2);
% end
% axis tight;

%% Save solution
% save('pvtol_ref_traj.mat','pos_coef','vel_coef','acc_coef','jer_coef','T_seg','coeff');

end

