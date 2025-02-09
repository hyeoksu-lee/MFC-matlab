function plot_jpdf_pres_vor(f_jpdf_pres_vor, pres, omega_xy)

    load variables/linestyle.mat;

    figure(f_jpdf_pres_vor);

    x = reshape(pres,[],1);
    y = reshape(omega_xy,[],1);

    [counts, xEdges, yEdges] = histcounts2(x, y, 100);

    % Convert histogram counts to probability density
    binWidthX = xEdges(2) - xEdges(1);
    binWidthY = yEdges(2) - yEdges(1);
    jointPDF = counts / (sum(counts(:)) * binWidthX * binWidthY);

    % Define bin centers
    xCenters = xEdges(1:end-1) + binWidthX/2;
    yCenters = yEdges(1:end-1) + binWidthY/2;

    % Plot joint PDF as a contour plot
    contourf(xCenters, yCenters, log(jointPDF'), 20, 'LineColor', 'none'); hold on;
    plot([pv pv],[0 5],'r--','LineWidth',1.5);
    xlim([-3 2]); xticks([-3:1:2]);
    ylim([0 5]); yticks([0:1:5]);
    colorbar; caxis([-10 6]);
    xlabel('$p$','Interpreter','latex');
    ylabel('$\omega_{xy}$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

end