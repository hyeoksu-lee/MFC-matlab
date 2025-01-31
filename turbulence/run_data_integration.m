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

    f1 = figure("DefaultAxesFontSize", 18);
    f2 = figure("DefaultAxesFontSize", 18);
    for l = 1:ncase
        load(post_stat_dir(l)+"/momentum_thickness.mat");

        mth_integrated = mth_integrated + mth / ncase;

        % For individual plot
        dmth = (mth(2:Nfiles(1)) - mth(1:Nfiles(1)-1)) ./ (time(2:Nfiles(1)) - time(1:Nfiles(1)-1));
        figure(f1); plot(time, mth, 'k-'); hold on; grid on;
        figure(f2); plot(time(1:Nfiles(1)-1), dmth, 'k-'); hold on; grid on;

    end

    % Integrated momentum thickness
    f_mom_thickness = figure("DefaultAxesFontSize", 18); 
    plot_mom_thickness(f_mom_thickness, time(1:Nfiles(1)), mth_integrated);
    saveas(f_mom_thickness, output_dir+"/momentum_thickness.png"); 
    close(f_mom_thickness); 
    disp("f_mom_thickness saved");

    % Integrated growth rate
    dmth = (mth_integrated(2:Nfiles(1)) - mth_integrated(1:Nfiles(1)-1)) ./ (time(2:Nfiles(1)) - time(1:Nfiles(1)-1));

    f_growth_rate = figure("DefaultAxesFontSize", 18); 
    plot_growth_rate(f_growth_rate, time(1:Nfiles(1)-1), dmth);
    saveas(f_growth_rate, output_dir+"/growth_rate.png"); 
    close(f_growth_rate); 
    disp("f_growth_rate saved");

    % Individual plot
    figure(f1);
    plot(time, mth_integrated, 'r-', 'LineWidth', 2); hold on; grid on;
    axis([0 500 0 14]);
    xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
    ylabel('$\delta / \delta_\omega^0$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    saveas(f1,output_dir+"/momentum_thickness_individual.png");
    close(f1);

    figure(f2);
    plot(time(1:Nfiles(1)-1), dmth, 'r-', 'LineWidth', 2); hold on; grid on;
    plot([time(1) time(end)],[0.014 0.014], 'b--', 'LineWidth', 2);
    axis([0 500 0 0.1]);
    xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
    ylabel('$\dot{\delta} / U_1$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    saveas(f2,output_dir+"/growth_rate_individual.png");
    close(f2);

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

    ruu_integrated = sqrt(ruu_integrated);
    rvv_integrated = sqrt(rvv_integrated);
    rww_integrated = sqrt(rww_integrated);
    ruv_integrated(ruv_integrated > 0) = 0;
    ruv_integrated = sqrt(-ruv_integrated);

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