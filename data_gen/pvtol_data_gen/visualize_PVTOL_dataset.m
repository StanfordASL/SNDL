%%
close all;

n_paths = 15; n_demos = 2;

cols = lines(15);
figure()
for i = 1:2:n_paths
    pi = p_idx_train(i);
    plot(paths_all_train{pi}(:,1),paths_all_train{pi}(:,2),'o-','Color',cols(i,:),...
     'markersize',15,'markerfacecolor',cols(i,:),'linewidth',2);
hold on
end
grid on
axis tight;
axis equal;
xlabel('p_x'); ylabel('p_z');
set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38)

[pos_coef,vel_coef,acc_coef,jer_coef,T_seg,coeff] = smooth_2D_path(paths{2},n_demos);
for i = 1:n_demos
    [x_d,y_d] = gen_poly(coeff{i},pos_coef,T_seg{i});
    plot(x_d,y_d,'m-','linewidth',2);
end
for i = 1:n_demos
    [Xs, Us, Ys] =  simulate_PVTOL(pos_coef,vel_coef,acc_coef,jer_coef,T_seg{i},coeff{i},0);
    plot(Xs(:,1),Xs(:,2),'r-','linewidth',2);
end
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%%
[pos_coef,vel_coef,acc_coef,jer_coef,T_seg,coeff] = smooth_2D_path(paths{2},n_demos);
for i = 1:n_demos
    [Xs, Us, Ys] =  simulate_PVTOL(pos_coef,vel_coef,acc_coef,jer_coef,T_seg{i},coeff{i},0);
    plot(Xs(:,1),Xs(:,2),'g-','linewidth',2);
end
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on
axis tight;
axis equal;

%%
[pos_coef,vel_coef,acc_coef,jer_coef,T_seg,coeff] = smooth_2D_path(paths{3},n_demos);
for i = 1:n_demos
    [Xs, Us, Ys] =  simulate_PVTOL(pos_coef,vel_coef,acc_coef,jer_coef,T_seg{i},coeff{i},0);
    plot(Xs(:,1),Xs(:,2),'r-','linewidth',2);
end
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on
axis tight;
axis equal;

 
    
function [x_d,y_d] = gen_poly(coeff,pos_coef,T_seg)

p_order = 10;
dT = 0.1;
break_T = [0;cumsum(T_seg)];
t_ref = (0:dT:break_T(end))';
x_d = zeros(length(t_ref),1); y_d = x_d; 

ns = 1;
for i = 1:length(t_ref)
    t = t_ref(i);
    if t > break_T(ns+1)
        ns = ns+1;
    end
    t_off = t - break_T(ns);
    pos_d = pos_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
    x_d(i) = pos_d(1); y_d(i) = pos_d(2);
    
end

end