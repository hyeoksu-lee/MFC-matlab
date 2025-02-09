function plot_min_pressure(f_mim_pressure, time, pres_min)

    load variables/linestyle.mat;
    figure(f_mim_pressure);

    plot(time,pres_min,'-ko','LineWidth',2); hold on; grid on;
    plot([0 400], [pv pv], 'r--','LineWidth',2);
    axis([0 400 -4 1]);
    yticks([-4:1:1]);
    xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
    ylabel('$p_{\mbox{min}}$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
end