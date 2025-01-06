function set_directory()

    % Number of cases
    ncase = 5;

    % Inputs
    mfc_dir = [ "/p/global/hyeoksu/MFC/reynolds/case/realization1/Re640/n511";
                "/p/global/hyeoksu/MFC/reynolds/case/realization2/Re640/n511";
                "/p/global/hyeoksu/MFC/reynolds/case/realization3/Re640/n511";
                "/p/global/hyeoksu/MFC/reynolds/case/realization4/Re640/n511";
                "/p/global/hyeoksu/MFC/reynolds/case/realization5/Re640/n511";
                ];

    % Outputs to input directory
    for i = 1:ncase
        if ~exist(strcat(mfc_dir(i),"/post_stat"), "dir")
            mkdir(strcat(mfc_dir(i),"/post_stat"));
        end
    end

    post_stat_dir = mfc_dir + "/post_stat";

    f_mean_streamwise_vel_dir = post_stat_dir + "/mean_streamwise_vel.png";
    f_Reynolds_stress_dir = post_stat_dir + "/Reynolds_stress.png";
    f_mom_thickness_dir = post_stat_dir + "/momentum_thickness.png";
    p_mom_thickness_dir = post_stat_dir + "/momentum_thickness.dat";

    f_bubble_radius_ratio_dir = post_stat_dir + "/bubble_radius_ratio_";
    f_bubble_radius_pdf_dir = post_stat_dir + "/bubble_radius_pdf_";
    p_bubble_radius_dir = post_stat_dir + "/bubble_radius_";

    % Outputs to this directory
    output_dir = "results/Re640/integrated";

    % Legend
    % legend_entity = ["$\tilde{M}_c = 0.1$", "$\tilde{M}_c = 0.05$", "$\tilde{M}_c = 0.01$", "$M_c = 0.0023$"];

    save variables/directory.mat;
end