function plot_growth_rate(f_growth_rate, time_in, dmth_in)

    load variables/linestyle.mat;
    figure(f_growth_rate);

    plot(time_in,dmth_in,'-ko'); hold on; grid on;
    plot([time_in(1) time_in(end)],[0.014 0.014], 'r--');
    axis([0 500 0 0.1]);
    xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
    ylabel('$\dot{\delta} / U_1$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
end