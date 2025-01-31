function plot_Reynolds_stress(f_Reynolds_stress, y_norm, ruu, rvv, rww, ruv)

    load variables/linestyle.mat;

    figure(f_Reynolds_stress);
    subplot(2,2,1); hold on;
    plot(y_norm(:,80:end),ruu(:,80:end),'k-'); grid on;
    axis([-8 8 0 0.3]);
    yticks([0:0.1:0.3]);
    xlabel('$y/\delta_m$','Interpreter','latex'); ylabel('$\sqrt{\left< \overline{\rho} u^{\prime 2} \right>} / U_1$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,2); hold on;
    plot(y_norm(:,80:end),rvv(:,80:end),'k-'); grid on;
    axis([-8 8 0 0.2]); 
    yticks([0:0.05:0.2]);
    xlabel('$y/\delta_m$','Interpreter','latex'); ylabel('$\sqrt{\left< \overline{\rho} v^{\prime 2} \right>} / U_1$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,3); hold on;
    plot(y_norm(:,80:end),rww(:,80:end),'k-'); grid on;
    axis([-8 8 0 0.2]); 
    yticks([0:0.05:0.2]);
    xlabel('$y/\delta_m$','Interpreter','latex'); ylabel('$\sqrt{\left< \overline{\rho} w^{\prime 2} \right>} / U_1$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,4); hold on;
    plot(y_norm(:,80:end),ruv(:,80:end),'k-'); grid on;
    axis([-8 8 0 0.2]); 
    yticks([0:0.05:0.2]);
    xlabel('$y/\delta_m$','Interpreter','latex'); ylabel('$\sqrt{-\left< \overline{\rho} u^{\prime} v^{\prime}\right>}  / U_1$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
end