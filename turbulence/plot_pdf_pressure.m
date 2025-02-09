function plot_pdf_pressure(f_pdf_pressure, pres)

    load variables/multiscale.mat;
    load variables/linestyle.mat;

    figure(f_pdf_pressure);
    histogram(reshape(pres,[],1),[-3:0.002:1.5],'EdgeColor','k','LineWidth',1.5,'Normalization','pdf', 'DisplayStyle', 'stairs'); hold on; grid on;
    plot([pv pv],[1e-6 1e4],'r--','LineWidth',1.5);
    xlim([-3 2]); 
    set(gca, 'YScale', 'log');
    xlabel('$p$','interpreter','latex');
    ylabel('$PDF$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

end