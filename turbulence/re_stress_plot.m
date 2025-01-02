% DIRECTORIES =============================================================================
ncase = 5;
mfc_loc = "../../delta/3d_mixing/apc_on/3d_mixing_phase1/";
% dirname = [ "Ma0_2/";"Ma0_005/";"Ma0_2_fine/";"Ma0_005_fine/";"Ma0_005/"];
dirname = [ "Ma0_2/";"Ma0_1/";"Ma0_05/";"Ma0_01/";"Ma0_005/"];
fresult = "./result/re_stress.png";
% =========================================================================================

colortype = ['b';'c';'g';'m';'r';'b';'g';'g';'r';'g';'r'];
linetype  = ["-";"-";"-";"-";"-";"-";"-";"-";"-";"-"];
markertype  = ["none";"none";"none";"none";"none";"^";"^";"^";"^"];
% markertype  = ["none";"none";"square";"square";"^";"^";"^";"^";"^"];

%% Import data
f=figure('DefaultAxesFontSize',18);
set(gcf,'Position',[100 100 1000 800]);
% A = readtable("./vreman97_ruu.dat");
% y = table2array(A(:,"Var1"));
% ruu = table2array(A(:,"Var2"));
% subplot(2,2,1);
% plot(y,ruu,'ko','LineWidth',1,'MarkerSize',4); hold on;

% A = readtable("./vreman97_rvv.dat");
% y = table2array(A(:,"Var1"));
% rvv = table2array(A(:,"Var2"));
% subplot(2,2,2);
% plot(y,rvv,'ko','LineWidth',1,'MarkerSize',4); hold on;

% A = readtable("./vreman97_rww.dat");
% y = table2array(A(:,"Var1"));
% rww = table2array(A(:,"Var2"));
% subplot(2,2,3);
% plot(y,rww,'ko','LineWidth',1,'MarkerSize',4); hold on;

% A = readtable("./vreman97_ruv.dat");
% y = table2array(A(:,"Var1"));
% ruv = table2array(A(:,"Var2"));
% subplot(2,2,4);
% plot(y,ruv,'ko','LineWidth',1,'MarkerSize',4); hold on;

for i=1:ncase
    fmthdata = mfc_loc + dirname(i) + "post_stat/re_stress_70.dat";
    A = readtable(fmthdata);
    y = table2array(A(:,"Var1"));
    ruu = table2array(A(:,"Var2"));
    rvv = table2array(A(:,"Var3"));
    rww = table2array(A(:,"Var4"));
    ruv = table2array(A(:,"Var5"));
    
    subplot(2,2,1);
    p1 = plot(y,ruu,'LineWidth',1,'MarkerSize',6); hold on;
    p1.Color = colortype(i); p1.LineStyle = linetype(i); p1.Marker= markertype(i);
    xlabel('$y$','interpreter','latex');
    ylabel('$\left< \overline{\rho} u^{\prime 2} \right>^{1/2}$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,2);
    p2 = plot(y,rvv,'LineWidth',1,'MarkerSize',6); hold on;
    p2.Color = colortype(i); p2.LineStyle = linetype(i); p2.Marker = markertype(i);
    xlabel('$y$','interpreter','latex');
    ylabel('$\left< \overline{\rho} v^{\prime 2} \right>^{1/2}$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,3);
    p3 = plot(y,rww,'LineWidth',1,'MarkerSize',6); hold on;
    p3.Color = colortype(i); p3.LineStyle = linetype(i); p3.Marker = markertype(i);
    xlabel('$y$','interpreter','latex');
    ylabel('$\left< \overline{\rho} w^{\prime 2} \right>^{1/2}$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,4);
    p4 = plot(y,ruv,'LineWidth',1,'MarkerSize',6); hold on;
    p4.Color = colortype(i); p4.LineStyle = linetype(i); p4.Marker = markertype(i);
    xlabel('$y$','interpreter','latex');
    ylabel('$-\left< \overline{\rho} u^{\prime} v^{\prime} \right>$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
end

% axis([0 100 0 8]); xticks([0:20:100]); yticks([0:2:8]);
% legend("Vreman et al. (1997), fDNS","Present, Ma = 0.2, 128^3","Present, Ma = 0.1, 128^3","Present, Ma = 0.05, 128^3","Present, Ma = 0.01, 128^3","Present, Ma = 0.005, 128^3",'location','eastoutside');
% legend("Vreman et al. (1996), 192^3","Vreman et al. (1996), 64^3","Present, 128^3",'location','northwest');
saveas(gcf,fresult);

disp("Done!");
exit;
