function plot_bubble_radius(f_bubble_radius_ratio, radius)

    load variables/poly_bubbles;
    
    figure(f_bubble_radius_ratio);
    semilogx(R0, radius./R0, 'k', 'LineWidth', 2); hold on;
    set(gca,'TickLabelInterpreter','latex');
    
end