function set_directory()

    % Number of cases
    ncase = 1;

    % Inputs
    mfc_dir = [ "/p/global/hyeoksu/MFC/reynolds/case/p009-ml3-004/Re160/sglrun";
                ];

    % Outputs to input directory
    for i = 1:ncase
        if ~exist(strcat(mfc_dir(i),"/post_stat"), "dir")
            mkdir(strcat(mfc_dir(i),"/post_stat"));
        end
        if ~exist(strcat(mfc_dir(i),"/post_stat/pdf_pressure"), "dir")
            mkdir(strcat(mfc_dir(i),"/post_stat/pdf_pressure"));
        end
        if ~exist(strcat(mfc_dir(i),"/post_stat/pdf_omega_xy"), "dir")
            mkdir(strcat(mfc_dir(i),"/post_stat/pdf_omega_xy"));
        end
        if ~exist(strcat(mfc_dir(i),"/post_stat/jpdf_pres_vor"), "dir")
            mkdir(strcat(mfc_dir(i),"/post_stat/jpdf_pres_vor"));
        end
    end

    post_stat_dir = mfc_dir + "/post_stat";

    f_mean_streamwise_vel_dir = post_stat_dir + "/mean_streamwise_vel.png";
    f_Reynolds_stress_dir = post_stat_dir + "/Reynolds_stress.png";
    f_mom_thickness_dir = post_stat_dir + "/momentum_thickness.png";
    p_mom_thickness_dir = post_stat_dir + "/momentum_thickness.dat";
    p_vor_thickness_dir = post_stat_dir + "/vorticity_thickness.dat";
    p_kolmogorov_length_dir = post_stat_dir + "/kolmogorov_length.dat";
    p_min_pressure_dir = post_stat_dir + "/min_pressure.dat";
    f_min_pressure_dir = post_stat_dir + "/min_pressure.png";
    f_pdf_pressure_dir = post_stat_dir + "/pdf_pressure/pdf_tstep_";
    f_pdf_omega_xy_dir = post_stat_dir + "/pdf_omega_xy/pdf_tstep_";
    f_jpdf_pres_vor_dir = post_stat_dir + "/jpdf_pres_vor/jpdf_tstep_";

    f_bubble_radius_ratio_dir = post_stat_dir + "/bubble_radius_ratio_";
    f_bubble_radius_pdf_dir = post_stat_dir + "/bubble_radius_pdf_";
    p_bubble_radius_dir = post_stat_dir + "/bubble_radius_";

    % Outputs to this directory
    output_dir = "results/p009-ml3-004/Re160/sglrun";

    save variables/directory.mat;
end