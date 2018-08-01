function analyze_PVTOL_results(N_tr,hybrid)

%N_tr: #points used in training (demonstration) set for regression loss
%hybrid: boolean to also present results for CCM-hybrid sims (if availabe)

%% Load results

load (strcat('results/PVTOL_UNCON_SIM_',num2str(N_tr),'.mat')); %Uncon dyn
load (strcat('results/PVTOL_UNCONU_SIM_',num2str(N_tr),'.mat')); %UnconU dyn

if (hybrid)
    load (strcat('results/PVTOL_CCM_SIM_',num2str(N_tr),'_1.mat')); %PVTOL_CCM_DATA (Hybrid)
    PVTOL_CCM_DATA_H = PVTOL_CCM_DATA; %change name
    clear PVTOL_CCM_DATA;
end
load (strcat('results/PVTOL_CCM_SIM_',num2str(N_tr),'_0.mat')); %PVTOL_CCM_DATA (LQR)

%t_ref,errors,t_span,X
N_sim = size(PVTOL_UNCON_DATA,1);

%% Uncon
err_norms_uncon = zeros(N_sim,1);
figure(1)
for i = 1:N_sim
    err_norms_uncon(i) = rms(norms(PVTOL_UNCON_DATA{i,2},2,2));
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
    err_norms_unconu(i) = rms(norms(PVTOL_UNCONU_DATA{i,2},2,2));
    plot(PVTOL_UNCONU_DATA{i,1}(1:end-1),norms(PVTOL_UNCONU_DATA{i,2},2,2),'linewidth',2);
    hold all;
end
title('UNCONU');
xlabel('Time [s]');
ylabel('$\|e\|$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% CCM_0
err_norms_CCM = zeros(N_sim,1);
figure(3)
for i = 1:N_sim
    err_norms_CCM(i) = rms(norms(PVTOL_CCM_DATA{i,2},2,2));
    plot(PVTOL_CCM_DATA{i,1}(1:end-1),norms(PVTOL_CCM_DATA{i,2},2,2),'linewidth',2);
    hold all;
end
title('CCM');
xlabel('Time [s]');
ylabel('$\|e\|$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% CCM_H
if (hybrid)
    err_norms_CCM_H = zeros(N_sim,1);
    figure(4)
    for i = 1:N_sim
        err_norms_CCM_H(i) = rms(norms(PVTOL_CCM_DATA_H{i,2},2,2));
        plot(PVTOL_CCM_DATA_H{i,1}(1:end-1),norms(PVTOL_CCM_DATA_H{i,2},2,2),'linewidth',2);
        hold all;
    end
    title('CCM - Hybrid');
    xlabel('Time [s]');
    ylabel('$\|e\|$','interpreter','latex');
    set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
end

%% Box-plot
figure()
plot(err_norms_unconu,'go','markerfacecolor','g','markersize',15);hold on
plot(err_norms_uncon,'bo','markerfacecolor','b','markersize',15); hold on
plot(err_norms_CCM,'ro','markerfacecolor','r','markersize',15);
legend('uncon (no reg)','uncon (l2 reg)', 'ccm reg');

figure()
h = boxplot([err_norms_uncon,err_norms_CCM],...
    'labels',{'R-R','CCM-R'},...
    'colors','br',...
    'OutlierSize',15);
set(h,'linew',2);
set(gca,'linew',2);
set(findall(gcf,'type','text'),'FontSize',38);set(gca,'FontSize',38);

if (hybrid)
    figure()
    plot(err_norms_CCM,'ro','markerfacecolor','r','markersize',15); hold on
    plot(err_norms_CCM_H,'mo','markerfacecolor','m','markersize',15);
    figure()
    boxplot([err_norms_CCM,err_norms_CCM_H],'labels',{'CCM','CCM-H'},'colors','rm');
end

%% traj comparisons - CCM
load(strcat('results/PVTOL_test_traj_ccm_',num2str(N_tr),'.mat'));

cols = lines(N_sim);
figure()
title('CCM');
hold on
for i = 1:N_sim
    plot(x_opt{i}(:,1),x_opt{i}(:,2),'color',cols(i,:),'linestyle','--','linewidth',2);
    plot(PVTOL_CCM_DATA{i,4}(:,1),PVTOL_CCM_DATA{i,4}(:,2),'color',cols(i,:),'linestyle','-','linewidth',2);
end
grid on

%vehicle overlayed plot
dT = 0.01/2;
dt = dT/2;
figure()
i = N_sim;
plot(x_opt{i}(:,1),x_opt{i}(:,2),'color','r','linestyle','--','linewidth',2);hold on
plot(PVTOL_CCM_DATA{i,4}(:,1),PVTOL_CCM_DATA{i,4}(:,2),'color','b','linestyle','-','linewidth',2);
for j = 1:(0.5/dT):size(x_opt{i},1)
        plot_pvtol(x_opt{i}(j,1),x_opt{i}(j,2),x_opt{i}(j,3),'r');
end
for j = 1:(0.5/dt):size(PVTOL_CCM_DATA{i,4},1)
    plot_pvtol(PVTOL_CCM_DATA{i,4}(j,1),PVTOL_CCM_DATA{i,4}(j,2),PVTOL_CCM_DATA{i,4}(j,3),'b');
end
plot(x_opt{i}(1,1),x_opt{i}(1,2),'go','markersize',10,'markerfacecolor','g')
plot(x_opt{i}(end,1),x_opt{i}(end,2),'ko','markersize',10,'markerfacecolor','k')
grid on
set(gca,'gridlinestyle','--')
title('CCM');

%% traj comparisons - CCM - H

if (hybrid)
    figure()
    title('CCM - Hybrid');
    hold on
    for i = 1:N_sim
        plot(x_opt{i}(:,1),x_opt{i}(:,2),'color',cols(i,:),'linestyle','--','linewidth',2);
        plot(PVTOL_CCM_DATA_H{i,4}(:,1),PVTOL_CCM_DATA_H{i,4}(:,2),'color',cols(i,:),'linestyle','-','linewidth',2);
    end
    grid on
end
clear x_opt;

%% traj comparisons - UNCON
load(strcat('results/PVTOL_test_traj_uncon_',num2str(N_tr),'.mat'));
figure()
title('Uncon');
hold on
for i = 1:N_sim
    plot(x_opt{i}(:,1),x_opt{i}(:,2),'color',cols(i,:),'linestyle','--','linewidth',2);
    plot(PVTOL_UNCON_DATA{i,4}(:,1),PVTOL_UNCON_DATA{i,4}(:,2),'color',cols(i,:),'linestyle','-','linewidth',2);
end
grid on

%vehicle overlayed plot
dT = 0.01/2;
dt = dT/2;
figure()
i = N_sim;
plot(x_opt{i}(:,1),x_opt{i}(:,2),'color','r','linestyle','--','linewidth',2);hold on
plot(PVTOL_UNCON_DATA{i,4}(:,1),PVTOL_UNCON_DATA{i,4}(:,2),'color','b','linestyle','-','linewidth',2);
for j = 1:(0.5/dT):size(x_opt{i},1)
        plot_pvtol(x_opt{i}(j,1),x_opt{i}(j,2),x_opt{i}(j,3),'r');
end
for j = 1:(0.5/dt):size(PVTOL_UNCON_DATA{i,4},1)
    plot_pvtol(PVTOL_UNCON_DATA{i,4}(j,1),PVTOL_UNCON_DATA{i,4}(j,2),PVTOL_UNCON_DATA{i,4}(j,3),'b');
end
plot(x_opt{i}(1,1),x_opt{i}(1,2),'go','markersize',10,'markerfacecolor','g')
plot(x_opt{i}(end,1),x_opt{i}(end,2),'ko','markersize',10,'markerfacecolor','k')
grid on
set(gca,'gridlinestyle','--')
title('Uncon');

clear x_opt;

%% traj comparisons - UNCONU
load(strcat('results/PVTOL_test_traj_unconu_',num2str(N_tr),'.mat'));
figure()
title('UnconU');
hold on
for i = 1:N_sim
    plot(x_opt{i}(:,1),x_opt{i}(:,2),'color',cols(i,:),'linestyle','--','linewidth',2);
    plot(PVTOL_UNCONU_DATA{i,4}(:,1),PVTOL_UNCONU_DATA{i,4}(:,2),'color',cols(i,:),'linestyle','-','linewidth',2);
end
grid on







