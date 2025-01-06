close all; clear all;
format long

% Remove pre-existing var files
if exist("variables", 'dir')
    delete variables/*.mat;
end

%% INPUTS
% Set directories
set_directory(); load variables/directory.mat;
% Set multiscale parameters
set_multiscale(ncase); load variables/multiscale.mat;
% Set grids
set_space(ncase); load variables/space.mat;
% Set time steps
set_time(ncase); load variables/time.mat;
% Set line style
set_linestyle(ncase); load variables/linestyle.mat;
% Set options
set_options(ncase); load variables/options.mat;
% Set polydisperse bubbles
set_poly_bubbles(); load variables/poly_bubbles.mat;
% Set index
set_index(); load variables/index.mat;

% Load data
%% POST-PROCESS

if (mom_thickness == "T")
    mth_integrated = zeros(1,Nfiles(1));

    for l = 1:ncase
        load(post_stat_dir(l)+"/momentum_thickness.mat");

        mth_integrated = mth_integrated + mth / ncase;

    end
    f_mom_thickness = figure("DefaultAxesFontSize", 18); 

    plot_mom_thickness(f_mom_thickness, time(1:Nfiles(1)), mth_integrated);

    saveas(f_mom_thickness, output_dir+"/momentum_thickness.png"); 
    close(f_mom_thickness); 
    disp("f_mom_thickness saved");
end

if (Reynolds_stress == "T")
    ruu_integrated = zeros(np(1),Nfiles(1));   % Reynolds stress: sqrt( rho u u)
    rvv_integrated = zeros(np(1),Nfiles(1));   % Reynolds stress: sqrt( rho v v)
    rww_integrated = zeros(np(1),Nfiles(1));   % Reynolds stress: sqrt( rho w w)
    ruv_integrated = zeros(np(1),Nfiles(1));   % Reynolds stress: sqrt(-rho u v)

    for l = 1:ncase

        load(post_stat_dir(l)+"/Reynolds_stress.mat");

        ruu_integrated = ruu_integrated + ruu / ncase;
        rvv_integrated = rvv_integrated + rvv / ncase;
        rww_integrated = rww_integrated + rww / ncase;
        ruv_integrated = ruv_integrated + ruv / ncase;

    end

    ruu_integrated = sqrt(ruu_integrated)/2;
    rvv_integrated = sqrt(rvv_integrated)/2;
    rww_integrated = sqrt(rww_integrated)/2;
    ruv_integrated = sqrt(-ruv_integrated)/2;

    f_Reynolds_stress = figure("DefaultAxesFontSize", 18); 
    plot_Reynolds_stress(f_Reynolds_stress, y_norm, ruu_integrated, rvv_integrated, rww_integrated, ruv_integrated);

    saveas(f_Reynolds_stress, output_dir+"/Reynolds_stress.png");
    close(f_Reynolds_stress);
    disp("f_Reynolds_stress saved");
end

% if (mean_streamwise_vel == "T")
%     u_mean_integrated = zeros(np(1),Nfiles(1));

%     for l = 1:ncase
%         load(post_stat_dir(l)+"/mean_streamwise_vel.mat");

%         u_mean_integrated = u_mean_integrated + u_mean / ncase;

%     end
%     f_mean_streamwise_vel = figure("DefaultAxesFontSize", 18); 

%     plot_mean_streamwise_vel(f_mean_streamwise_vel, u_mean_integrated, y_norm);

%     saveas(f_mean_streamwise_vel, output_dir+"/mean_streamwise_vel.png"); 
%     close(f_mean_streamwise_vel); 
%     disp("f_mean_streamwise_vel saved");
% end