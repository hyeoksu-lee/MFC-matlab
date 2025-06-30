
function compare_mom_thickness()

    load variables/directory.mat;
    load variables/time.mat;
    load variables/linestyle.mat;

    
    f_compare_mom_thickness = figure("DefaultAxesFontSize", 18);
    set(f_compare_mom_thickness,"position",[100 100 700 600]);

    A = readmatrix("data/vreman96.dat");
    t = A(:,1);
    mth = A(:,2);
    plot(t, mth, 'ko', 'LineWidth', 2); hold on; grid on;

    for l = 1:ncase
        A = readmatrix(p_mom_thickness_dir(l));
        t = A(:,1);
        mth = A(:,2);
        plot(t, mth, linestyle(l), 'LineWidth', 2); hold on; grid on;
        ylim([0.5 0.9]);
        yticks([0.5:0.1:0.9]);
    end
    legend("Vreman et al. (1996), Ma = 0.2","one-way","two-way",'interpreter','latex','location','northwest');
    xlabel("$t$",'interpreter','latex');
    ylabel("$\delta_m$",'interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    xlim([0 15]);
    saveas(f_compare_mom_thickness,output_dir+"/f_compare_mom_thickness.png"); 
    close(f_compare_mom_thickness);
    disp('compare_mom_thickness saved');

end