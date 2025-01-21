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

%% POST-PROCESS
parpool('local');
for l = 1:ncase

    % Preallocate arrays
    mth = zeros(1,Nfiles(l));       % Momentum thickness
    ruu = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho u u)
    rvv = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho v v)
    rww = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho w w)
    ruv = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho u v)
    u_mean = zeros(np(l),Nfiles(l));% Mean streamwise velocity
    y_norm = zeros(np(l),Nfiles(l));

    % Parallized loop over timesteps
    parfor q = 1:Nfiles(l)
        % Read data
        filename = strcat(mfc_dir(l),"/restart_data/lustre_",int2str(timesteps(l,q)),".dat");
        disp(filename);
        fileID = fopen(filename,'r');
        A = fread(fileID,'double');
        fclose(fileID);

        % Reassign density & velocity components
        rho = reshape(A(1:r(l)), mp(l), np(l), pp(l)); % Continuity
        mom = reshape(A(r(l)+1:4*r(l)), 3, mp(l), np(l), pp(l)); % Momentum
        vel = mom ./ reshape(rho, 1, mp(l), np(l), pp(l)); % Velocity

        clear A

        % Compute mean and fluctuating quantities
        rho_mean = mean(rho,[1 3]); % Mean density
        mom_mean = squeeze(mean(mom,[2 4])); % Mean momentum
        vel_mean = mom_mean./rho_mean; % Favre-averaged velocity

        % Compute momentum thickness
        f = rho_mean .* (1 - mom_mean(1, :) ./ rho_mean) .* (1 + mom_mean(1, :) ./ rho_mean) / 4;
        mth(q) = dy(l)*trapz(f);

        % Normalized y-axis 
        y_norm(:,q) = y(l,1:np(l))'/mth(q);

        % Mean streamwise velocity
        u_mean(:,q) = vel_mean(1,:)';

        if (Reynolds_stress == "T")
            % Compute Reynolds stress
            vel_fluc = vel - permute(repmat(vel_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]);
            ruu(:, q) = mean(rho .* vel_fluc(1, :, :, :).^2, [1 3]);
            rvv(:, q) = mean(rho .* vel_fluc(2, :, :, :).^2, [1 3]);
            rww(:, q) = mean(rho .* vel_fluc(3, :, :, :).^2, [1 3]);
            ruv(:, q) = mean(rho .* vel_fluc(1, :, :, :) .* vel_fluc(2, :, :, :), [1 3]);
        end
    end

    if (mean_streamwise_vel == "T")
        save_mean_streamwise_vel(post_stat_dir(l), time(l,1:Nfiles(l)), y_norm, u_mean); disp('mean_streamwise_vel data saved');
    end
    if (Reynolds_stress == "T")
        save_Reynolds_stress(post_stat_dir(l), time(l,1:Nfiles(l)), y_norm, ruu, rvv, rww, ruv); disp('Reynolds_stress data saved');
    end
    if (mom_thickness == "T")
        f_mom_thickness = figure("DefaultAxesFontSize", 18); plot_mom_thickness(f_mom_thickness, time(l,1:Nfiles(l)), mth);
        saveas(f_mom_thickness, f_mom_thickness_dir(l)); close(f_mom_thickness); disp('f_mom_thickness saved');
        save_mom_thickness(post_stat_dir(l), time(l,1:Nfiles(l)), mth); disp('mom_thickness data saved');
        print_mom_thickness(p_mom_thickness_dir(l), time(l,1:Nfiles(l)), mth); disp('p_momentum_thickness printed');
    end
end

% compare_mom_thickness();
% compare_bubble_radius_pdf();
% compare_bubble_radius_ratio();