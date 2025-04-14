function set_directory()

    %% Number of cases
    ncase = 3;

    %% Proj dir
    proj_dir = ["weno3m/n127";
                "weno5m/n127";
                "wcns6ld/n127";
                ];

    %% Input data location
    mfc_dir = "/p/global/hyeoksu/MFC/reynolds/case/3d_taylorgreen/"+proj_dir;

    %% Outputs to matlab directory
    output_dir = "results/n127/compare.png";

    %% Create post_stat directory
    for i = 1:ncase
        if ~exist(strcat(mfc_dir(i),"/post_stat"), "dir")
            mkdir(strcat(mfc_dir(i),"/post_stat"));
        end
    end

    save variables/directory.mat;
end