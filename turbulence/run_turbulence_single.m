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
for l = 1:ncase

    % Preallocate arrays
    vth = zeros(1,Nfiles(l));       % Vorticity thickness
    mth = zeros(1,Nfiles(l));       % Momentum thickness
    ruu = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho u u)
    rvv = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho v v)
    rww = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho w w)
    ruv = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho u v)
    u_mean = zeros(np(l),Nfiles(l));% Mean streamwise velocity
    y_norm = zeros(np(l),Nfiles(l));
    eta_min = zeros(1,Nfiles(l));
    eta_max = zeros(1,Nfiles(l));

    % Parallized loop over timesteps
    for q = 1:Nfiles(l)
        % Read data
        filename = strcat(mfc_dir(l),"/restart_data/lustre_",int2str(timesteps(l,q)),".dat");
        disp(filename);
        fileID = fopen(filename,'r');
        A = fread(fileID,'double');
        fclose(fileID);

        % Reassign density & velocity components
        rho = reshape(A(1:r(l)), mp(l), np(l), pp(l)); % Continuity
        mom = permute(reshape(A(r(l)+1:4*r(l)), mp(l), np(l), pp(l), 3), [4 1 2 3]) * (u_b(l)/u_m(l)); % Momentum
        vel = mom ./ repmat(reshape(rho, 1, mp(l), np(l), pp(l)), [3 1 1 1]); % Velocity
        
        E = reshape(A((E_idx - 1)*r(l) + 1:E_idx*r(l)), mp(l), np(l), pp(l)) * (u_b(l)/u_m(l))^2; % Energy
        dyn_p = 0.5 * squeeze(mom(1,:,:,:).*vel(1,:,:,:) + mom(2,:,:,:).*vel(2,:,:,:) + mom(3,:,:,:).*vel(3,:,:,:));
        pres = (1/Gamma_a) * (E - dyn_p) - Pi_inf_a(l);

        clear A, E, dyn_p

        % Compute mean and fluctuating quantities
        rho_mean = mean(rho,[1 3]); % Mean density
        mom_mean = squeeze(mean(mom,[2 4])); % Mean momentum
        vel_mean = mom_mean./rho_mean; % Favre-averaged velocity
        pres_mean = mean(pres, [1 3]);

        % Compute vorticity thickness
        dudy = (vel_mean(2:end) - vel_mean(1:end-1)) / dy(l);
        dudymax = max(abs(dudy));
        vth(q) = 2 / dudymax;

        % Compute momentum thickness
        f = rho_mean .* (1 - vel_mean(1, :)) .* (1 + vel_mean(1, :)) / 4;
        mth(q) = dy(l)*trapz(f);

        % Normalized y-axis 
        y_norm(:,q) = y(l,1:np(l))'/mth(q);

        % Mean streamwise velocity
        u_mean(:,q) = vel_mean(1,:)';


        if (Reynolds_stress == "T")
            % Compute velocity and pressure fluctuation
            vel_fluc = vel - permute(repmat(vel_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]);

            % Compute Reynolds stress
            ruu(:, q) = mean(rho .* squeeze(vel_fluc(1, :, :, :)).^2, [1 3]);
            rvv(:, q) = mean(rho .* squeeze(vel_fluc(2, :, :, :)).^2, [1 3]);
            rww(:, q) = mean(rho .* squeeze(vel_fluc(3, :, :, :)).^2, [1 3]);
            ruv(:, q) = mean(rho .* squeeze(vel_fluc(1, :, :, :)) .* squeeze(vel_fluc(2, :, :, :)), [1 3]);

            dudx = squeeze((vel_fluc(1, 3:end, 2:end-1, 2:end-1) - vel_fluc(1, 1:end-2, 2:end-1, 2:end-1)) / (2*dx(l)));
            dvdx = squeeze((vel_fluc(2, 3:end, 2:end-1, 2:end-1) - vel_fluc(2, 1:end-2, 2:end-1, 2:end-1)) / (2*dx(l)));
            dwdx = squeeze((vel_fluc(3, 3:end, 2:end-1, 2:end-1) - vel_fluc(3, 1:end-2, 2:end-1, 2:end-1)) / (2*dx(l)));

            dudy = squeeze((vel_fluc(1, 2:end-1, 3:end, 2:end-1) - vel_fluc(1, 2:end-1, 1:end-2, 2:end-1)) / (2*dy(l)));
            dvdy = squeeze((vel_fluc(2, 2:end-1, 3:end, 2:end-1) - vel_fluc(2, 2:end-1, 1:end-2, 2:end-1)) / (2*dy(l)));
            dwdy = squeeze((vel_fluc(3, 2:end-1, 3:end, 2:end-1) - vel_fluc(3, 2:end-1, 1:end-2, 2:end-1)) / (2*dy(l)));

            dudz = squeeze((vel_fluc(1, 2:end-1, 2:end-1, 3:end) - vel_fluc(1, 2:end-1, 2:end-1, 1:end-2)) / (2*dy(l)));
            dvdz = squeeze((vel_fluc(2, 2:end-1, 2:end-1, 3:end) - vel_fluc(2, 2:end-1, 2:end-1, 1:end-2)) / (2*dy(l)));
            dwdz = squeeze((vel_fluc(3, 2:end-1, 2:end-1, 3:end) - vel_fluc(3, 2:end-1, 2:end-1, 1:end-2)) / (2*dy(l)));

            % Strain-rate tensor
            s11 = dudx;     
            s12 = 0.5*(dudy + dvdx);    
            s13 = 0.5*(dudz + dwdx);

            s21 = s12;
            s22 = dvdy;
            s23 = 0.5*(dvdz + dwdy);

            s31 = s13;
            s32 = s23;
            s33 = dwdz;

            ss = s11.^2 + s12.^2 + s13.^2 + s21.^2 + s22.^2 + s23.^2 + s31.^2 + s32.^2 + s33.^2;

            % Dissipation rate
            eps_mean = (2/Re0)*mean(ss, [1 3]);
            
            % Kolmogorov length scale
            eta_max(q) = (1/(Re0^3 * min(eps_mean)))^(0.25);
            eta_min(q) = (1/(Re0^3 * max(eps_mean)))^(0.25);



            % Pressure statistics
            pres_fluc = pres - permute(repmat(pres_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]);

        end
    end

    if (mean_streamwise_vel == "T")
        save_mean_streamwise_vel(post_stat_dir(l), time(l,1:Nfiles(l)), y_norm, u_mean); disp('mean_streamwise_vel data saved');
    end
    if (Reynolds_stress == "T")
        save_Reynolds_stress(post_stat_dir(l), time(l,1:Nfiles(l)), y_norm, ruu, rvv, rww, ruv); disp('Reynolds_stress data saved');
        print_kolmogorov_length(p_kolmogorov_length_dir(l), time(l,1:Nfiles(l)), eta_max, eta_min); disp('p_kolmogorov_length printed');
    end
    if (mom_thickness == "T")
        f_mom_thickness = figure("DefaultAxesFontSize", 18); plot_mom_thickness(f_mom_thickness, time(l,1:Nfiles(l)), mth);
        saveas(f_mom_thickness, f_mom_thickness_dir(l)); close(f_mom_thickness); disp('f_mom_thickness saved');
        save_mom_thickness(post_stat_dir(l), time(l,1:Nfiles(l)), mth); disp('mom_thickness data saved');
        print_mom_thickness(p_mom_thickness_dir(l), time(l,1:Nfiles(l)), mth); disp('p_momentum_thickness printed');

        % Vorticity thickness
        print_vor_thickness(p_vor_thickness_dir(l), time(l,1:Nfiles(l)), vth); disp('p_vorticity_thickness printed');
    end
end

% compare_mom_thickness();
% compare_bubble_radius_pdf();
% compare_bubble_radius_ratio();