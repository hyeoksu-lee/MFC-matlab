function plot_bubble_radius_pdf(f_bubble_radius_pdf, radius)

    load variables/poly_bubbles;
    
    figure(f_bubble_radius_pdf);
    loglog(radius, weight, 'k', 'LineWidth', 2); hold on;
    loglog(R0, weight, 'k:', 'LineWidth', 1);
    set(gca,'TickLabelInterpreter','latex');

end