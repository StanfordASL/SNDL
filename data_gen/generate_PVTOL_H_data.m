clear; close all; clc


%% Circle data


X_all = [];
U_all = [];
Y_all = [];

load('pvtol_h_data/circle_cw.mat');

X_all = [X_all; X];
U_all = [U_all; U];
Y_all = [Y_all; Y];

clear X U Y;

load('pvtol_h_data/circle_cw_fast.mat');

X_all = [X_all; X];
U_all = [U_all; U];
Y_all = [Y_all; Y];

clear X U Y;

load('pvtol_h_data/circle_cw_2fast.mat');

X_all = [X_all; X];
U_all = [U_all; U];
Y_all = [Y_all; Y];

clear X U Y;

load('pvtol_h_data/circle_ccw.mat');

X_all = [X_all; X];
U_all = [U_all; U];
Y_all = [Y_all; Y];

clear X U Y;

load('pvtol_h_data/circle_ccw_fast.mat');

X_all = [X_all; X];
U_all = [U_all; U];
Y_all = [Y_all; Y];

clear X U Y;

load('pvtol_h_data/circle_ccw_2fast.mat');

X_all = [X_all; X];
U_all = [U_all; U];
Y_all = [Y_all; Y];

clear X U Y;

%% Swoop data

X_swoop = [];
U_swoop = [];
Y_swoop = [];

load('pvtol_h_data/swoop_fwd.mat');

X_swoop = [X_swoop; X];
U_swoop = [U_swoop; U];
Y_swoop = [Y_swoop; Y];

clear X U Y;

load('pvtol_h_data/swoop_bwd.mat');

X_swoop = [X_swoop; X];
U_swoop = [U_swoop; U];
Y_swoop = [Y_swoop; Y];

clear X U Y;

%% Plot data

figure()
subplot(5,1,1)
histogram(X_all(:,2));
title('z');

subplot(5,1,2)
histogram(X_all(:,3));
title('\phi');

subplot(5,1,3)
histogram(X_all(:,4));
title('vy');

subplot(5,1,4)
histogram(X_all(:,5));
title('vz')

subplot(5,1,5)
histogram(X_all(:,6));
title('\phi dot');

%% Generate a random sample too

N_halton = 500;

raw_range = [-5,-3,-25*pi/180,-2,-2,-0.8;
              5,0 , 25*pi/180, 2, 2, 0.8];

rand_halton = generateHaltonSamples(6,N_halton);
rand_halton = repmat(raw_range(1,:),N_halton,1) + rand_halton.*repmat(raw_range(2,:)-raw_range(1,:),N_halton,1);

%% Plot random samples

figure()
subplot(5,1,1)
histogram(rand_halton(:,2));
title('z');

subplot(5,1,2)
histogram(rand_halton(:,3));
title('\phi');

subplot(5,1,3)
histogram(rand_halton(:,4));
title('vy');

subplot(5,1,4)
histogram(rand_halton(:,5));
title('vz')

subplot(5,1,5)
histogram(rand_halton(:,6));
title('\phi dot');

%% Separate into train, val, and test

ind_all = (1:size(X_all,1))';

ind_train = ind_all( floor(linspace(1,size(X_all,1),2000) ) ); 

ind_rem = setdiff(ind_all,ind_train);

ind_val = ind_rem( floor(linspace(1,size(ind_rem,1),3000) ) );

ind_rem = setdiff(ind_rem, ind_val);

ind_test = ind_rem( floor(linspace(1,size(ind_rem,1),2000)) );


ind_all_swoop = (1:size(X_swoop,1))';

ind_train_swoop = ind_all_swoop( floor(linspace(1,size(X_swoop,1),600) ) ); 

ind_rem_swoop = setdiff(ind_all_swoop,ind_train_swoop);

ind_val_swoop = ind_rem_swoop( floor(linspace(1,size(ind_rem_swoop,1),700) ) );

ind_rem_swoop = setdiff(ind_rem_swoop, ind_val_swoop);

ind_test_swoop = ind_rem_swoop( floor(linspace(1,size(ind_rem_swoop,1),700)) );

%% Assemble Train, Val, and Test sets

N_max = length(ind_train) + length(ind_train_swoop);

X = [X_all(ind_train,:);
     X_swoop(ind_train_swoop,:);
     rand_halton];
 
U = [U_all(ind_train,:);
     U_swoop(ind_train_swoop,:)];

Y = [Y_all(ind_train,:);
     Y_swoop(ind_train_swoop,:)];

save('PVTOL_H_data_train.mat','X','U','Y','N_max');

X_val = [X_all(ind_val,:);
         X_swoop(ind_val_swoop,:)];

U_val = [U_all(ind_val,:);
         U_swoop(ind_val_swoop,:)];

Y_val = [Y_all(ind_val,:);
         Y_swoop(ind_val_swoop,:)];

save('PVTOL_H_data_val.mat','X_val','U_val','Y_val');

X_test = [X_all(ind_test,:);
         X_swoop(ind_test_swoop,:)];

U_test = [U_all(ind_test,:);
         U_swoop(ind_test_swoop,:)];

Y_test = [Y_all(ind_test,:);
         Y_swoop(ind_test_swoop,:)];

save('PVTOL_H_data_test.mat','X_test','U_test','Y_test');










