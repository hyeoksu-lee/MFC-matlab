function set_directory()

    % Number of cases
    ncase = 1;

    % Proj dir
    proj_dir = ["p010-ml3-001/Re50/weno_Re_flux_F/weno-cu6/N511";
                ];

    % Inputs
    mfc_dir = "/p/global/hyeoksu/MFC/reynolds/case/"+proj_dir;

    % Outputs to matlab directory
    % output_dir = "results/p010-ml3-001/Re50/weno_Re_flux_T";
    output_dir = "results/"+proj_dir;

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
        if ~exist(strcat(mfc_dir(i),"/post_stat/jpdf_pres_omega_xy"), "dir")
            mkdir(strcat(mfc_dir(i),"/post_stat/jpdf_pres_omega_xy"));
        end
        if ~exist(strcat(mfc_dir(i),"/post_stat/tke_budget"), "dir")
            mkdir(strcat(mfc_dir(i),"/post_stat/tke_budget"));
        end
        if ~exist(strcat(mfc_dir(i),"/post_stat/energy_spectrum"), "dir")
            mkdir(strcat(mfc_dir(i),"/post_stat/energy_spectrum"));
        end
    end

    save variables/directory.mat;
end