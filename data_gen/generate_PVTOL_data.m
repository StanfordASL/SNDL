function generate_PVTOL_data()
clc; close all;

addpath('pvtol_data_gen');

raw_range = [7,7,75*pi/180,5,5,3*pi/4];

%% Pre-gen set of fixed waypoint paths

N_paths = 30;
paths_all_train = create_2D_path(N_paths);

%treat above paths as centers of Gaussian blobs

%% Train data

n_traj = 30; %#trajectories recorded
n_demos = 2; %#demonstrations per unique path
n_p = n_traj/n_demos;

[X, U, Y, p_idx_train] = get_dataset_pvtol(paths_all_train,n_p,n_demos);
N_max = size(X,1); %max # of (x,u,x_dot) tuples

%additinally add some Halton sampled points
rand_halton = generateHaltonSamples(6,500);
rand_halton = repmat(-raw_range,500,1) + 2*rand_halton.*repmat(raw_range,500,1);
X = [X;rand_halton]; 

save('PVTOL_data_train.mat','X','U','Y','paths_all_train','p_idx_train','N_max');

%% Val data

paths_all_val = create_2D_path(N_paths);

n_traj = 30; %#trajectories recorded
n_demos = 2; %#demonstrations per unique path
n_p = n_traj/n_demos;

[X_val, U_val, Y_val, p_idx_val] = get_dataset_pvtol(paths_all_val,n_p,n_demos);

save('PVTOL_data_val.mat','X_val','U_val','Y_val','paths_all_val','p_idx_val');

%% Final plot

figure(3)
plot(X(:,1),X(:,2),'bx');hold on
plot(X_val(:,1),X_val(:,2),'ro');
xlabel('x'); ylabel('z');

figure(4)
plot3(X(:,3),X(:,4),X(:,5),'bx');hold on
plot3(X_val(:,3),X_val(:,4),X_val(:,5),'ro');
xlabel('\phi'); ylabel('vx'); zlabel('vz');


end

function [X, U, Y, p_idx] = get_dataset_pvtol(paths_all,n_p,n_demos)

wp_var = 0.1;

X = [];
U = [];
Y = [];

%first randomly select waypoint paths
p_idx = randperm(length(paths_all),n_p);

n_seg = size(paths_all{1},1);

for i = 1:n_p
    %randomly distort path
    wp_path = paths_all{p_idx(i)} + mvnrnd([0;0],wp_var*eye(2),n_seg);
        
    %smooth and simulate
    for j = 1:n_demos
        sim_ok = 0;
        while(~sim_ok)
            fprintf('(%d,%d)\n',i,j);
            %smooth
            [pos_coef,vel_coef,acc_coef,jer_coef,T_seg,coeff] = smooth_2D_path(wp_path,1);
            %simulate
            [Xs, Us, Ys] =  simulate_PVTOL(pos_coef,vel_coef,acc_coef,jer_coef,T_seg{1},coeff{1});
            %check sanity
            if max(norms(Xs(:,1:2),2,2))<=10 && max(abs(Xs(:,3)))<=80*(pi/180)
                sim_ok = 1;
            end
        end
        X = [X; Xs];
        U = [U; Us];
        Y = [Y; Ys];
    end
end

end