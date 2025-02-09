function plot_mom_thickness(f_mom_thickness, time_in, mth_in)

    load variables/linestyle.mat;
    figure(f_mom_thickness);

    plot(time_in,mth_in,'-ko','LineWidth',2); hold on; grid on;
    axis([0 400 0 14]);
    xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
    ylabel('$\delta / \delta_\omega^0$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
end