
function compare_bubble_radius_pdf()

    load variables/directory.mat;
    load variables/time.mat;
    load variables/linestyle.mat;

    for q = 1:length(Nt_compare)
        f_compare_bubble_radius_pdf = figure("DefaultAxesFontSize", 18);
        set(f_compare_bubble_radius_pdf,"position",[100 100 600 600]);
        for l = 1:ncase
            A = readmatrix(strcat(p_bubble_radius_dir(l),int2str(Nt_compare(q)),".dat"));
            R0 = A(:,1);
            radius = A(:,2);
            weight = A(:,3);
            loglog(radius, weight, linestyle(l), 'LineWidth', 2); hold on; grid on;
        end
        loglog(R0, weight, 'k:', 'LineWidth', 1);

        xlabel("$R/R_{0ref}$",'interpreter','latex');
        ylabel("$f(R)$",'interpreter','latex');
        time_on_title = "$t = "+num2str(time(1,q))+"$";
        title(time_on_title,'interpreter','latex');
        set(gca,'TickLabelInterpreter','latex');
        saveas(f_compare_bubble_radius_pdf,output_dir+"/compare_bubble_radius_pdf/"+int2str(Nt_compare(q))+".png"); 
        close(f_compare_bubble_radius_pdf);
        disp('compare_bubble_radius_pdf saved');
    end

end