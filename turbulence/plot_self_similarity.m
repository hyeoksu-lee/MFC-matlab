function plot_self_similarity(f_self_similarity, u_mean, y_norm, nt)

    load variables/linestyle.mat;

    figure(f_self_similarity);
    plot(u_mean,y_norm,'color',blueGrad(nt,:)); hold on;
    % axis([-1.1 1.1 -5 5]);
    xlabel('$\bar{u}/\delta$','Interpreter','latex');
    ylabel('y/\delta');
    % xticks([-1 -0.5 0 0.5 1]);
    % yticks([-5:5]);
    set(gca,'TickLabelInterpreter','latex');

end