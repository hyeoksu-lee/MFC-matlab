function plot_mean_streamwise_vel(f_mean_streamwise_vel, u_mean, y_norm)

    load variables/linestyle.mat;

    figure(f_mean_streamwise_vel);
    plot(u_mean,y_norm,'k-'); hold on;
    xlim([-1.1 1.1]);
    ylim([-10 10]);
    xlabel('$\bar{u}/\delta$','Interpreter','latex');
    ylabel('y/\delta');
    % xticks([-1 -0.5 0 0.5 1]);
    % yticks([-5:5]);
    set(gca,'TickLabelInterpreter','latex');

end