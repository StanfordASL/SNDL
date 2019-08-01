function analyze_PVTOL_results(N_tr,ctrl)
close all;

%N_tr: #points used in training (demonstration) set for regression loss
%ctrl: controller type {lqr, mpc, ccm}

%% Load results

folder = 'results';
load (strcat(folder,'_new','/PVTOL_CCM_',num2str(N_tr),'_',upper(ctrl),'_new','.mat')); %CCM dyn

%use lqr controller for other models
if strcmp(ctrl,'ccm')
    ctrl = 'lqr';
end
load (strcat(folder,'/PVTOL_UNCON_',num2str(N_tr),'_',upper(ctrl),'.mat')); %Uncon dyn
load (strcat(folder,'/PVTOL_UNCONU_',num2str(N_tr),'_',upper(ctrl),'.mat')); %UnconU dyn

%t_ref,errors,t_span,X
N_sim = size(PVTOL_UNCON_DATA,1);

%% Uncon
err_norms_uncon = zeros(N_sim,1);
figure(1)
for i = 1:N_sim
%     err_norms_uncon(i) = rms(norms(PVTOL_UNCON_DATA{i,2},2,2));
    t_ref_i = PVTOL_UNCON_DATA{i,1}(1:end-1);
    err_norm_i = norms(PVTOL_UNCON_DATA{i,2},2,2);
    err_norms_uncon(i) = sqrt( (1/t_ref_i(end)) * trapz( t_ref_i, err_norm_i.^2));
    plot(PVTOL_UNCON_DATA{i,1}(1:end-1),norms(PVTOL_UNCON_DATA{i,2},2,2),'linewidth',2);
    hold all;
end
title('UNCON');
xlabel('Time [s]');
ylabel('$\|e\|$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% UnconU
err_norms_unconu = zeros(N_sim,1);
figure(2)
for i = 1:N_sim
%     err_norms_unconu(i) = rms(norms(PVTOL_UNCONU_DATA{i,2},2,2));
    t_ref_i = PVTOL_UNCONU_DATA{i,1}(1:end-1);
    err_norm_i = norms(PVTOL_UNCONU_DATA{i,2},2,2);
    err_norms_unconu(i) = sqrt( (1/t_ref_i(end)) * trapz( t_ref_i, err_norm_i.^2));
    plot(PVTOL_UNCONU_DATA{i,1}(1:end-1),norms(PVTOL_UNCONU_DATA{i,2},2,2),'linewidth',2);
    hold all;
end
title('UNCONU');
xlabel('Time [s]');
ylabel('$\|e\|$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% CCM
err_norms_CCM = zeros(N_sim,1);
figure(3)
for i = 1:N_sim
%     err_norms_CCM(i) = rms(norms(PVTOL_CCM_DATA{i,2},2,2));
    t_ref_i = PVTOL_CCM_DATA{i,1}(1:end-1);
    err_norm_i = norms(PVTOL_CCM_DATA{i,2},2,2);
    err_norms_CCM(i) = sqrt( (1/t_ref_i(end)) * trapz( t_ref_i, err_norm_i.^2));
    plot(PVTOL_CCM_DATA{i,1}(1:end-1),norms(PVTOL_CCM_DATA{i,2},2,2),'linewidth',2);
    hold all;
end
title('CCM');
xlabel('Time [s]');
ylabel('$\|e\|$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Box-plot
figure(4)
plot(err_norms_unconu,'go','markerfacecolor','g','markersize',15);hold on
plot(err_norms_uncon,'bo','markerfacecolor','b','markersize',15); hold on
plot(err_norms_CCM,'ro','markerfacecolor','r','markersize',15);
legend('uncon (no reg)','uncon (l2 reg)', 'ccm reg');
ylabel('$\|e\|$','interpreter','latex');

figure(5)
h = boxplot([err_norms_unconu,err_norms_uncon,err_norms_CCM],...
    'labels',{'N-R','R-R','CCM-R'},...
    'colors','mbr',...
    'OutlierSize',15);
set(h,'linew',2);
set(gca,'linew',2);
set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38);

keyboard;

%% traj comparisons - CCM
load(strcat(folder,'_new','/PVTOL_test_traj_ccm_',num2str(N_tr),'.mat'));
cols = lines(N_sim);
figure(10)
title('CCM');
hold on
for i = 1:5:N_sim
    plot(x_opt{i}(:,1),x_opt{i}(:,2),'color',cols(i,:),'linestyle','--','linewidth',2);
    plot(PVTOL_CCM_DATA{i,4}(:,1),PVTOL_CCM_DATA{i,4}(:,2),'color',cols(i,:),'linestyle','-','linewidth',2);
end
grid on

keyboard;


% %% Find worst trajs for UNCON
% 
% [uncon_worst, worst_idx] = sort(err_norms_CCM,'descend');
% worst_idx = worst_idx(1:5)
% 
% %% traj comparisons - CCM
% clear x_opt;
% 
% load(strcat(folder,'/PVTOL_test_traj_ccm_',num2str(N_tr),'.mat'));
% 
% cols = lines(length(worst_idx));
% figure(6)
% title('CCM');
% hold on
% for j = 1:length(worst_idx)
%     i = worst_idx(j);
%     plot(x_opt{i}(:,1),x_opt{i}(:,2),'color',cols(j,:),'linestyle','--','linewidth',1.5);
%     plot(PVTOL_CCM_DATA{i,4}(:,1),PVTOL_CCM_DATA{i,4}(:,2),'color',cols(j,:),'linestyle','-','linewidth',1.5);
%     scatter(x_opt{i}(1,1),x_opt{i}(1,2),300,cols(j,:),'filled');
%     scatter(PVTOL_CCM_DATA{i,4}(end,1),PVTOL_CCM_DATA{i,4}(end,2),300,'r','filled');
% end
% grid on
% % axis equal
% axis tight
% set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38);
% 
% %vehicle overlayed plot
% dT = 0.01/2;
% dt = dT/2;
% figure(7)
% subplot(2,1,1)
% i = worst_idx(1);
% plot(x_opt{i}(:,1),x_opt{i}(:,2),'color','r','linestyle','--','linewidth',2);hold on
% plot(PVTOL_CCM_DATA{i,4}(:,1),PVTOL_CCM_DATA{i,4}(:,2),'color','b','linestyle','-','linewidth',2);
% for j = 1:(0.5/dT):size(x_opt{i},1)
%         plot_pvtol(x_opt{i}(j,1),x_opt{i}(j,2),x_opt{i}(j,3),'r');
% end
% for j = 1:(0.5/dt):size(PVTOL_CCM_DATA{i,4},1)
%     plot_pvtol(PVTOL_CCM_DATA{i,4}(j,1),PVTOL_CCM_DATA{i,4}(j,2),PVTOL_CCM_DATA{i,4}(j,3),'b');
% end
% plot(x_opt{i}(1,1),x_opt{i}(1,2),'go','markersize',10,'markerfacecolor','g')
% plot(x_opt{i}(end,1),x_opt{i}(end,2),'ko','markersize',10,'markerfacecolor','k')
% plot(PVTOL_CCM_DATA{i,4}(end,1),PVTOL_CCM_DATA{i,4}(end,2),'bo','markersize',10,'markerfacecolor','b')
% grid on
% set(gca,'gridlinestyle','--')
% title('CCM');
% set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38);
% % axis equal 
% axis tight
% 
% clear x_opt;
% 
% %% traj comparisons - UNCON
% load(strcat(folder,'/PVTOL_test_traj_uncon_',num2str(N_tr),'.mat'));
% figure(8)
% title('Uncon');
% hold on
% for j = 1:length(worst_idx)
%     i = worst_idx(j);
%     plot(x_opt{i}(:,1),x_opt{i}(:,2),'color',cols(j,:),'linestyle','--','linewidth',1.5);
%     plot(PVTOL_UNCON_DATA{i,4}(:,1),PVTOL_UNCON_DATA{i,4}(:,2),'color',cols(j,:),'linestyle','-','linewidth',1.5);
%     scatter(x_opt{i}(1,1),x_opt{i}(1,2),300,cols(j,:),'filled');
%     scatter(PVTOL_UNCON_DATA{i,4}(end,1),PVTOL_UNCON_DATA{i,4}(end,2),300,'r','filled');
% end
% grid on
% % axis equal
% axis tight
% set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38);
% 
% %vehicle overlayed plot
% dT = 0.01/2;
% dt = dT/2;
% figure(7)
% subplot(2,1,2)
% i = worst_idx(1);
% plot(x_opt{i}(:,1),x_opt{i}(:,2),'color','r','linestyle','--','linewidth',2);hold on
% plot(PVTOL_UNCON_DATA{i,4}(:,1),PVTOL_UNCON_DATA{i,4}(:,2),'color','b','linestyle','-','linewidth',2);
% for j = 1:(0.5/dT):size(x_opt{i},1)
%         plot_pvtol(x_opt{i}(j,1),x_opt{i}(j,2),x_opt{i}(j,3),'r');
% end
% for j = 1:(0.5/dt):size(PVTOL_UNCON_DATA{i,4},1)
%     plot_pvtol(PVTOL_UNCON_DATA{i,4}(j,1),PVTOL_UNCON_DATA{i,4}(j,2),PVTOL_UNCON_DATA{i,4}(j,3),'b');
% end
% plot(x_opt{i}(1,1),x_opt{i}(1,2),'go','markersize',10,'markerfacecolor','g')
% plot(x_opt{i}(end,1),x_opt{i}(end,2),'ko','markersize',10,'markerfacecolor','k')
% plot(PVTOL_UNCON_DATA{i,4}(end,1),PVTOL_UNCON_DATA{i,4}(end,2),'bo','markersize',10,'markerfacecolor','b')
% grid on
% set(gca,'gridlinestyle','--')
% title('Uncon');
% set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38);
% % axis equal
% axis tight
% 
% clear x_opt;

% %% traj comparisons - UNCONU
% load(strcat(folder,'/PVTOL_test_traj_unconu_',num2str(N_tr),'.mat'));
% cols = lines(N_sim);
% figure(10)
% title('UnconU');
% hold on
% for i = 1:5:N_sim
%     plot(x_opt{i}(:,1),x_opt{i}(:,2),'color',cols(i,:),'linestyle','--','linewidth',2);
%     plot(PVTOL_UNCONU_DATA{i,4}(:,1),PVTOL_UNCONU_DATA{i,4}(:,2),'color',cols(i,:),'linestyle','-','linewidth',2);
% end
% grid on







