% DIRECTORIES =============================================================================
ncase = 4;
mfc_loc = "../../delta/3d_mixing/apc_on/3d_mixing_phase1/";
dirname = [ "Ma0_2/";"Ma0_005/";"Ma0_2_fine/";"Ma0_005_fine/";"Ma0_005/"];
fresult = "./result/mth.png";
% =========================================================================================

colortype = ['b';'r';'b';'r';'r';'b';'g';'g';'r';'g';'r'];
linetype  = ["-";"-";":";":";"-";"-";"-";"-";"-";"-"];
markertype  = ["none";"none";"square";"square";"-";"-";"-";"-";"-";"-"];

%% Import data
f=figure('DefaultAxesFontSize',18);
% set(gcf,'Position',[100 100 350 300]);
set(gcf,'Position',[100 100 300 450]);

mth = zeros(ncase,21);
for i=1:ncase
    fmthdata = mfc_loc + dirname(i) + "post_stat/mom_thickness.dat";
    A = readtable(fmthdata);
    time = table2array(A(:,"Var1"));
    mth(i,:) = table2array(A(:,"Var2"));
    p = plot(time,mth(i,:),'LineWidth',1,'MarkerSize',4.5); hold on;
    p.Color = colortype(i);
    p.LineStyle = linetype(i);
    p.Marker = markertype(i);
end

% Vreman's data ===========================================================================
A = readtable("./vreman96_grid192_mth.dat");
time = table2array(A(:,"Var1"));
mth = table2array(A(:,"Var2"));
plot(time,mth,'ko','LineWidth',1,'MarkerSize',4.5); hold on;

A = readtable("./vreman96_grid64_mth.dat");
time = table2array(A(:,"Var1"));
mth = table2array(A(:,"Var2"));
plot(time,mth,'k^','LineWidth',1,'MarkerSize',4.5); hold on;
% =========================================================================================

xlabel('$tU_1/ \delta^{0}_{\omega}$','interpreter','latex');
ylabel('$\delta/ \delta^{0}_{\omega}$','interpreter','latex');
axis([0 80 0 6]); xticks([0:20:80]); yticks([0:1:6]); 
legend( "$\mbox{Vreman et al. (1996), Ma = 0.2, }192^3$", ... 
        "$\mbox{Vreman et al. (1996), Ma = 0.2, }64^3$", ...
        "$\mbox{Present, Ma = 0.2, }128^3$", ...
        "$\mbox{Present, Ma = 0.005, }128^3$", ...
        "$\mbox{Present, Ma = 0.2, }192^3$", ...
        "$\mbox{Present, Ma = 0.005, }192^3$", ...
        'location','southoutside','interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
saveas(gcf,fresult);

disp("Done!");
exit;
