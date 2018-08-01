function [Xs, Us, Ys] =  simulate_PVTOL(pos_coef,vel_coef,acc_coef,jer_coef,T_seg,coeff)

%% load trajectory

% load('pvtol_ref_traj.mat');

p_order = 10;
dT = 0.01;

%Create splines for each segment
break_T = [0;cumsum(T_seg)];
t_ref = (0:dT:break_T(end))';
x_d = zeros(length(t_ref),1); y_d = x_d; 
vx_d = x_d; vy_d = x_d;
ax_d = x_d; ay_d = x_d;
jx_d = x_d; jy_d = x_d;
ns = 1;
for i = 1:length(t_ref)
    t = t_ref(i);
    if t > break_T(ns+1)
        ns = ns+1;
    end
    t_off = t - break_T(ns);
    pos_d = pos_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
    x_d(i) = pos_d(1); y_d(i) = pos_d(2);
    
    vel_d = vel_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
    vx_d(i) = vel_d(1); vy_d(i) = vel_d(2);
    
    acc_d = acc_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
    ax_d(i) = acc_d(1); ay_d(i) = acc_d(2);
    
    jer_d = jer_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
    jx_d(i) = jer_d(1); jy_d(i) = jer_d(2);
end

X_ref = [x_d,y_d,vx_d,vy_d,ax_d,ay_d,jx_d];

%% Simulate

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
X0 = [x_d(1),y_d(1),0,0,0,0];
dt = 0.01;
t_end = t_ref(end);
t_span = (0:dt:t_end)';
X = zeros(length(t_span),6);
U = zeros(length(t_span)-1,2);
Y = zeros(length(t_span)-1,6);

X(1,:) = X0;

fprintf('Starting sim...');
for i = 1:length(t_span)-1
    %get ctrl
    Xd = interp1(t_ref,X_ref,t_span(i));
    u = pvtol_ctrl(Xd,X(i,:)');
    U(i,:) = u';
    %get dx_dt for data
    dX_dt = pvtol_dyn(0,X(i,:)',u);
    Y(i,:) = dX_dt';
    
    %simulate
%     X(i+1,:) = X(i,:) + Y(i,:)*dt;
    [~, xs] = ode113(@(ts,xs) pvtol_dyn(ts,xs,u),...
                     [t_span(i),t_span(i+1)],X(i,:)',ode_options);
                 
    X(i+1,:) = xs(end,:);  
    X(i+1,3) = wrapToPi(X(i+1,3));
end
fprintf('Done!\n');

Xs = X(1:10:end-1,:);
Us = U(1:10:end,:);
Ys = Y(1:10:end,:);

%% Record and save

% figure(3)
% title('reference and simulated paths');
% hold on
% plot(X_ref(:,1),X_ref(:,2),'r-','linewidth',2);
% plot(X(:,1),X(:,2),'b-','linewidth',2);
% 
% keyboard;

end

function u = pvtol_ctrl(Xd,X)
mass = 0.486;
g = 9.81;
len = 0.25;
kp = 2; kd = 0.1;

x = X(1); y = X(2); phi = wrapToPi(X(3));
vx = X(4); vy = X(5); phi_dot = X(6);

x_dot = vx*cos(phi) - vy*sin(phi);
y_dot = vx*sin(phi) + vy*cos(phi);

u_thrust = mass*(g+Xd(6) + kd*(Xd(4)-y_dot) + kp*(Xd(2)-y));

x_ddot = -(u_thrust/mass)*sin(phi);

phi_c = (-1/g)*(Xd(5) + kd*(Xd(3)-x_dot) + kp*(Xd(1)-x));
phi_c_dot = wrapToPi((-1/g)*(Xd(7) + kd*(Xd(5)-x_ddot) + kp*(Xd(3)-x_dot)));

u_torque = kp*(phi_c-phi) + kd*(phi_c_dot-phi_dot);

u = [1, 1;
    len,-len]\[u_thrust;u_torque];

%clip control
for i = 1:2
    u(i) = min(max(u(i),1.5),10);
end

end


