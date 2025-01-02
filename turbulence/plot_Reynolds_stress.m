function plot_Reynolds_stress(f_Reynolds_stress, y_norm, ruu, rvv, rww, ruv, nt)

    load variables/linestyle.mat;

    figure(f_Reynolds_stress);
    subplot(2,2,1); hold on;
    plot(y_norm,ruu,'color',blueGrad(nt,:));
    % axis([-5 5 0 0.5]); 
    xlabel('$y/\delta_m$','Interpreter','latex'); ylabel('$\left< \overline{\rho} u^{\prime 2} \right>^{1/2}$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,2); hold on;
    plot(y_norm,rvv,'color',blueGrad(nt,:));
    % axis([-5 5 0 0.5]); 
    xlabel('$y/\delta_m$','Interpreter','latex'); ylabel('$\left< \overline{\rho} v^{\prime 2} \right>^{1/2}$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,3); hold on;
    plot(y_norm,rww,'color',blueGrad(nt,:));
    % axis([-5 5 0 0.5]); 
    xlabel('$y/\delta_m$','Interpreter','latex'); ylabel('$\left< \overline{\rho} w^{\prime 2} \right>^{1/2}$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    subplot(2,2,4); hold on;
    plot(y_norm,ruv,'color',blueGrad(nt,:));
    % axis([-5 5 -0.2 0]); 
    xlabel('$y/\delta_m$','Interpreter','latex'); ylabel('$-\left< \overline{\rho} u^{\prime} v^{\prime}\right>$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
end