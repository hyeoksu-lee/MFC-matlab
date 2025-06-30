function plot_growth_rate(f_growth_rate, time_in, dmth_in)

    load variables/linestyle.mat;
    figure(f_growth_rate);

    plot(time_in,dmth_in,'-ko','LineWidth',2); hold on; grid on;
    plot([0 400],[0.014 0.014], 'r--','LineWidth',2);
    plot([0 400],[0.016 0.016], 'b--','LineWidth',2);
    plot([0 400],[0.022 0.022], 'g--','LineWidth',2);
    axis([0 400 0 0.1]);
    xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
    ylabel('$\dot{\delta} / U_1$','interpreter','latex');
    legend("$\mbox{Present}$","$0.014^{[1-4]}$","$0.016^{[5-6]}$","$0.022^{[1]}$",'Interpreter','latex','location','northeast');
    set(gca,'TickLabelInterpreter','latex');
    
end