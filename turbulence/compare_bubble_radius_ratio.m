
function compare_bubble_radius_ratio()

    load variables/directory.mat;
    load variables/time.mat;
    load variables/linestyle.mat;

    for q = 1:length(Nt_compare)
        f_compare_bubble_radius_ratio = figure("DefaultAxesFontSize", 18);
        set(f_compare_bubble_radius_ratio,"position",[100 100 600 600]);
        for l = 1:ncase
            A = readmatrix(strcat(p_bubble_radius_dir(l),int2str(Nt_compare(q)),".dat"));
            R0 = A(:,1);
            radius = A(:,2);
            weight = A(:,3);
            loglog(R0, radius./R0, linestyle(l), 'LineWidth', 2); hold on; grid on;
            ylim([1 1e2]);
        end
        legend("one-way","two-way",'interpreter','latex','location','northwest');
        xlabel("$R_{0}$",'interpreter','latex');
        ylabel("$R/R_{0}$",'interpreter','latex');
        time_on_title = "$t = "+num2str(time(1,q))+"$";
        title(time_on_title,'interpreter','latex');
        set(gca,'TickLabelInterpreter','latex');
        saveas(f_compare_bubble_radius_ratio,output_dir+"/compare_bubble_radius_ratio/"+int2str(Nt_compare(q))+".png"); 
        close(f_compare_bubble_radius_ratio);
        disp('compare_bubble_radius_ratio saved');
    end

end