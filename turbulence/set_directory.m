function set_directory()

    % Number of cases
    ncase = 2;

    % Inputs
    mfc_dir = [ "/p/work1/hyeoksu/MFC/one-way-aM/case/aM/two-way";
                "/p/work1/hyeoksu/MFC/one-way-aM/case/aM/one-way";
                ];

    % Outputs to input directory
    % for i = 1:ncase
    %     if ~exist(strcat(mfc_dir(i),"/post_stat"), "dir")
    %         mkdir(strcat(mfc_dir(i),"/post_stat"));
    %     end
    % end
    f_self_similarity_dir = mfc_dir + "/post_stat/self_similarity.png";
    f_Reynolds_stress_dir = mfc_dir + "/post_stat/Reynolds_stress.png";
    f_mom_thickness_dir = mfc_dir + "/post_stat/momentum_thickness.png";
    p_mom_thickness_dir = mfc_dir + "/post_stat/momentum_thickness.dat";
    f_bubble_radius_ratio_dir = mfc_dir + "/post_stat/bubble_radius_ratio_";
    f_bubble_radius_pdf_dir = mfc_dir + "/post_stat/bubble_radius_pdf_";
    p_bubble_radius_dir = mfc_dir + "/post_stat/bubble_radius_";

    % Outputs to this directory
    output_dir = "./results/aps";

    % Legend
    % legend_entity = ["$\tilde{M}_c = 0.1$", "$\tilde{M}_c = 0.05$", "$\tilde{M}_c = 0.01$", "$M_c = 0.0023$"];

    save variables/directory.mat;
end