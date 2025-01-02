function plot_mom_thickness(f_mom_thickness, time_in, mth_in)

    load variables/linestyle.mat;
    figure(f_mom_thickness);

    plot(time_in,mth_in,'-ko'); hold on;
    % axis([0 250 0 8]);
    xlabel('t');
    ylabel('$\delta$','interpreter','latex');
    % xticks([0:50:250]);
    % yticks([0:1:8]);
    set(gca,'TickLabelInterpreter','latex');
    
end