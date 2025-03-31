function plot_mean_streamwise_vel(f_mean_streamwise_vel, u_mean, y_norm)

    load variables/linestyle.mat;

    figure(f_mean_streamwise_vel);
    A = readmatrix("data/Wang_et_al_2022/umean.dat");
    y_norm_ref1 = A(:,1); umean_ref1 = A(:,2);
    plot(umean_ref1,y_norm_ref1,'ro','LineWidth',2); hold on;
    plot(u_mean(:,25:end)/2,y_norm(:,25:end),'k-');
    xlim([-0.5 0.5]);
    ylim([-5 5]);
    xlabel('$\bar{u}/\Delta U$','Interpreter','latex');
    ylabel('y/\delta');
    legend("Wang et al. (2022)", "Present",'interpreter','latex','location','northwest');
    set(gca,'TickLabelInterpreter','latex');

end