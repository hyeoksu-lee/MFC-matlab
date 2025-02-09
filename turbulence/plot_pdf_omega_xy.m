function plot_pdf_omega_xy(f_pdf_omega_xy, omega_xy)

    load variables/linestyle.mat;

    figure(f_pdf_omega_xy);
    histogram(reshape(omega_xy,[],1),'EdgeColor','k','LineWidth',1.5,'Normalization','pdf', 'DisplayStyle', 'stairs'); hold on; grid on;
    xlim([0 5]); 
    set(gca, 'YScale', 'log');
    xlabel('$\omega_{xy}$','interpreter','latex');
    ylabel('$PDF$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

end