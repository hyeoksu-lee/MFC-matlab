close all; clear all;
format long

% Remove pre-existing var files
if exist("variables", 'dir')
    delete variables/*.mat;
end

%% INPUTS
disp("initialize")
tic
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
toc

%% POST-PROCESS
for l = 1:ncase
    disp("setup each case")
    tic
    set_global_parameters(l); load variables/global_parameters.mat;

    set_grid(); load variables/grid.mat;

    set_allocate_temporal_data(); load variables/temporal_data.mat;
    toc

    % Parallized loop over timesteps
    % parpool("Processes",2)
    if (read_raw_data == "T")
        for q = 1:Nfiles
            
            % task = getCurrentTask();
            % fprintf('Worker %d is processing file index %d\n', task.ID, q);

            % Read data in conservative form
            disp("read data")
            tic
            qc = f_read_data(strcat(mfc_dir(l),"/restart_data/lustre_",int2str(timesteps(q)),".dat"));
            toc

            % Convert conservative variables to primitive variables
            disp("convert_cons_to_prim");
            tic
            qp = f_convert_cons_to_prim(qc);
            toc

            % Compute mean quantities
            disp("qp_mean_fluc");
            tic
            [qp_mean qp_fluc] = f_compute_qp_mean_fluc(qp);
            toc

            % Save mean streamwise velocity
            u_mean(:,q) = qp_mean(momxb,:)';

            if (pres_stat == "T")
                disp("pres_stat")
                tic
                % Compute minimum pressure PDF of minimum pressure
                pres_min(q) = min(qp(E_idx,:,:,:),[],"all");
                % Compute PDF of pressure
                plot_pdf_pressure(qp(E_idx,:,:,:), timesteps(q));
                toc
            end

            % Compute derivatives
            disp("vel_derivatives")
            tic
            [dvel_ds dvelmean_dy dvelfluc_ds] = f_compute_vel_derivatives(qp(momxb:momxe,:,:,:), qp_mean(momxb:momxe,:), qp_fluc(momxb:momxe,:,:,:));
            dpmean_dy = f_compute_derivative_1d(squeeze(qp_mean(E_idx,:)),y_cc);
            toc

            if (vorticity == "T")
                disp("vorticity")
                tic
                % Compute vorticity
                omega = f_compute_vorticity(dvel_ds);
                omega_xy = sqrt(omega(1,:,:,:).^2 + omega(2,:,:,:).^2);
                plot_pdf_omega_xy(omega_xy,timesteps(q));
                if (pres_stat == "T")
                    plot_jpdf_pres_omega_xy(qp(E_idx,:,:,:),omega_xy,timesteps(q));
                end
                toc
            end

            % Compute mixing layer thickness
            disp("mixing_layer_thickness")
            tic
            [vth(q) mth(q) y_norm_mth(:,q) y_norm_vth(:,q)] = f_compute_mixing_layer_thickness(qp_mean,dvelmean_dy);
            toc

            if (Reynolds_stress == "T")
                % Compute Reynolds stress
                disp("Reynolds_stress")
                tic
                [ruu(:,q) rvv(:,q) rww(:,q) ruv(:,q)] = f_compute_Reynolds_stress(qp(1,:,:,:), qp_fluc(momxb:momxe,:,:,:));
                plot_Reynolds_stress(y_norm_vth(:,q), ruu(:,q), rvv(:,q), rww(:,q), ruv(:,q), timesteps(q));
                toc
            end

            if (kolmogorov == "T")
                % Compute Kolmogorov length scale / dy_min
                disp("kolmogorov")
                tic
                eta_min(q) = f_compute_kolmogorov_scale(dvelfluc_ds);
                toc
            end

            if (tke_budget == "T")
                % Compute TKE budget
                disp("TKE budget")
                tic
                [T P D S C] = f_compute_tke_budget(dvel_ds, dvelmean_dy, dvelfluc_ds, dpmean_dy, squeeze(qp(1,:,:,:)), qp_fluc(momxb:momxe,:,:,:), squeeze(qp(E_idx,:,:,:)), mth(q), y_norm_mth(:,q), timesteps(q));
                plot_tke_budget(y_norm_mth(:,q), T, P, D, S, C, timesteps(q));
                toc
                % Compute integrated TKE budget
                tic
                Pintegrated(q) = trapz(y_cc,P);
                Dintegrated(q) = trapz(y_cc,D);
                toc
            end

            if (energy_spectrum == "T")
                % Compute 1D energy spectrum
                disp("energy spectrum")
                tic
                [k,E] = f_compute_energy_spectrum(qp_fluc, mth(q));
                plot_energy_spectrum(k, E, timesteps(q));
                toc
            end

        end

        disp("save, plot, print")
        tic
        % SAVE DATA AS .MAT FILES
        save(post_stat_dir+"/mean_streamwise_vel.mat","time","y_norm_mth","u_mean");
        save(post_stat_dir+"/mixlayer_thickness.mat","time","mth","vth");
        if (Reynolds_stress == "T")
            save(post_stat_dir+"/Reynolds_stress.mat","time","y_norm_mth","ruu","rvv","rww","ruv");
        end
 
        % SAVE DATA AS .DAT FILES
        if (pres_stat == "T")
            print_min_pressure(time, pres_min);
        end
        print_mixlayer_thickness(time, mth, vth);
        if (kolmogorov == "T")
            print_kolmogorov_length(time, eta_min);
        end

        % PLOT
        plot_mom_thickness(time, mth);
        plot_growth_rate(time, mth);
        if (pres_stat == "T")
            plot_min_pressure(time, pres_min);
        end
        if (Reynolds_stress == "T")
            plot_Reynolds_stress(y_norm_vth, ruu, rvv, rww, ruv, "all");
        end
        if (tke_budget == "T")
            plot_integrated_tke_budget(time, Pintegrated, Dintegrated);
        end
        toc
    end

    if (tke_budget == "T" && avg_over_self_sim == "T")
        disp("average_tke_budget")
        tic
        f_average_tke_budget();
        toc
    end
end

%% FUNCTIONS
% f_read_data
function qc = f_read_data(filename)

    load variables/global_parameters.mat;
    load variables/index.mat;

    % Read data
    disp(filename);
    fileID = fopen(filename,'r');
    A = fread(fileID,'double');
    fclose(fileID);

    % Reassign density & velocity components
    qc = permute(reshape(A, mp, np, pp, sys_size),[4 1 2 3]);
    qc(momxb:momxe,:,:,:) = qc(momxb:momxe,:,:,:) * (u_b/u_m);  % Momentum
    qc(E_idx,:,:,:) = qc(E_idx,:,:,:) * (u_b/u_m)^2;    % Energy

end

% f_convert_cons_to_prim
function qp = f_convert_cons_to_prim(qc)

    load variables/index.mat;
    load variables/multiscale.mat;

    qp = zeros(size(qc));

    % density
    qp(1,:,:,:) = qc(1,:,:,:);

    dyn_p = 0;
    for i = momxb:momxe
        qp(i,:,:,:) = qc(i,:,:,:) ./ qc(1,:,:,:); % Velocity
        dyn_p = dyn_p + 0.5*qc(i,:,:,:).*qp(i,:,:,:); % Dynamic pressure
    end

    % pressure
    qp(E_idx,:,:,:) = (1/Gamma)*(qc(E_idx,:,:,:) - dyn_p - Pi_inf);
end

% f_compute_qp_mean_fluc
function [qp_mean qp_fluc] = f_compute_qp_mean_fluc(qp)

    load variables/global_parameters.mat;
    load variables/index.mat

    % Compute mean quantities
    qp_mean = zeros(sys_size,np);

    % mean rho
    qp_mean(1,:) = squeeze(mean(qp(1,:,:,:),[2 4]));

    % favre-averaged velocity
    for i = momxb:momxe
        qp_mean(i,:) = squeeze(mean(qp(1,:,:,:).*qp(i,:,:,:),[2 4])./mean(qp(1,:,:,:),[2 4]));
    end

    % mean pressure
    qp_mean(E_idx,:) = squeeze(mean(qp(E_idx,:,:,:),[2 4]));

    % density fluctuation
    qp_fluc = qp - permute(repmat(qp_mean, [1 1 mp pp]), [1 3 2 4]);
end

% f_compute_vel_derivatives
function [dvel_ds dvelmean_dy dvelfluc_ds] = f_compute_vel_derivatives(vel, vel_mean, vel_fluc)

    load variables/global_parameters.mat;
    load variables/grid.mat;
    load variables/options.mat
    
    dvel_ds = zeros(3,3,mp,np,pp);
    dvelmean_dy = zeros(3,np);
    dvelfluc_ds = zeros(3,3,mp,np,pp);

    % Compute velocity derivatives
    if (vorticity == "T" || tke_budget == "T")
        vel1 = squeeze(vel(1,:,:,:));
        vel2 = squeeze(vel(2,:,:,:));
        vel3 = squeeze(vel(3,:,:,:));
    
        dvel_ds(1,1,:,:,:) = f_compute_derivative_3d(vel1,x_cc,1);
        dvel_ds(2,1,:,:,:) = f_compute_derivative_3d(vel2,x_cc,1);
        dvel_ds(3,1,:,:,:) = f_compute_derivative_3d(vel3,x_cc,1);

        dvel_ds(1,2,:,:,:) = f_compute_derivative_3d(vel1,y_cc,2);
        dvel_ds(2,2,:,:,:) = f_compute_derivative_3d(vel2,y_cc,2);
        dvel_ds(3,2,:,:,:) = f_compute_derivative_3d(vel3,y_cc,2);

        dvel_ds(1,3,:,:,:) = f_compute_derivative_3d(vel1,z_cc,3);
        dvel_ds(2,3,:,:,:) = f_compute_derivative_3d(vel2,z_cc,3);
        dvel_ds(3,3,:,:,:) = f_compute_derivative_3d(vel3,z_cc,3);
    end

    % favre-averaged velocity derivatives
    dvelmean_dy(1,:) = f_compute_derivative_1d(vel_mean(1,:),y_cc);
    dvelmean_dy(2,:) = f_compute_derivative_1d(vel_mean(2,:),y_cc);
    dvelmean_dy(3,:) = f_compute_derivative_1d(vel_mean(3,:),y_cc);

    % fluctuating velocity derivatives
    if (kolmogorov == "T" || tke_budget == "T")
        vel_fluc1 = squeeze(vel_fluc(1,:,:,:));
        vel_fluc2 = squeeze(vel_fluc(2,:,:,:));
        vel_fluc3 = squeeze(vel_fluc(3,:,:,:));
        
        dvelfluc_ds(1,1,:,:,:) = f_compute_derivative_3d(vel_fluc1,x_cc,1);
        dvelfluc_ds(2,1,:,:,:) = f_compute_derivative_3d(vel_fluc2,x_cc,1);
        dvelfluc_ds(3,1,:,:,:) = f_compute_derivative_3d(vel_fluc3,x_cc,1);

        dvelfluc_ds(1,2,:,:,:) = f_compute_derivative_3d(vel_fluc1,y_cc,2);
        dvelfluc_ds(2,2,:,:,:) = f_compute_derivative_3d(vel_fluc2,y_cc,2);
        dvelfluc_ds(3,2,:,:,:) = f_compute_derivative_3d(vel_fluc3,y_cc,2);

        dvelfluc_ds(1,3,:,:,:) = f_compute_derivative_3d(vel_fluc1,z_cc,3);
        dvelfluc_ds(2,3,:,:,:) = f_compute_derivative_3d(vel_fluc2,z_cc,3);
        dvelfluc_ds(3,3,:,:,:) = f_compute_derivative_3d(vel_fluc3,z_cc,3);
    end
end

% f_compute_mixing_layer_thickness
function [vth mth y_norm_mth y_norm_vth] = f_compute_mixing_layer_thickness(qp_mean, dvelmean_dy)

    load variables/grid.mat;
    load variables/index.mat;

    % Compute vorticity thickness
    vth = 2 / max(abs(dvelmean_dy(1,:)),[],"all");
    
    % Compute momentum thickness
    f = qp_mean(1,:) .* (1 - qp_mean(momxb, :)) .* (1 + qp_mean(momxb, :)) / 4;
    mth = trapz(y_cc,f);

    % Compute mth-normalized y value
    y_norm_mth = y_cc'/mth;
    y_norm_vth = y_cc'/vth;

end

% f_compute_Reynolds_stress
function [ruu rvv rww ruv] = f_compute_Reynolds_stress(rho, vel_fluc)

    rho = squeeze(rho);
    vel_fluc1 = squeeze(vel_fluc(1, :, :, :));
    vel_fluc2 = squeeze(vel_fluc(2, :, :, :));
    vel_fluc3 = squeeze(vel_fluc(3, :, :, :));

    % Compute Reynolds stress
    ruu = mean(rho .* vel_fluc1.^2, [1 3]) / 2;
    rvv = mean(rho .* vel_fluc2.^2, [1 3]) / 2;
    rww = mean(rho .* vel_fluc3.^2, [1 3]) / 2;
    ruv = mean(rho .* vel_fluc1.*vel_fluc2, [1 3]) / 2;
end

% f_compute_strain_rate_tensor
function [SS S11 S12 S13 S22 S23 S33] = f_compute_strain_rate_tensor(dvel_ds);
    S11 = squeeze(dvel_ds(1,1,:,:,:));
    S12 = squeeze(0.5*(dvel_ds(1,2,:,:,:) + dvel_ds(2,1,:,:,:)));
    S13 = squeeze(0.5*(dvel_ds(1,3,:,:,:) + dvel_ds(3,1,:,:,:)));

    S22 = squeeze(dvel_ds(2,2,:,:,:));
    S23 = squeeze(0.5*(dvel_ds(2,3,:,:,:) + dvel_ds(3,2,:,:,:)));

    S33 = squeeze(dvel_ds(3,3,:,:,:));

    SS = S11.^2 + S12.^2 + S13.^2 + S12.^2 + S22.^2 + S23.^2 + S13.^2 + S23.^2 + S33.^2;
end

% f_compute_kolmogorov_scale
function eta_min = f_compute_kolmogorov_scale(dvelfluc_ds)

    load variables/global_parameters.mat;
    load variables/grid.mat;

    % Compute turbulent strain-rate tensor
    [ss s11 s12 s13 s21 s22 s33] = f_compute_strain_rate_tensor(dvelfluc_ds);

    % dissipation
    eps_mean = (2/Re)*mean(ss, [1 3]);
    eta_min = (1/(Re^3 * max(eps_mean)))^(0.25) /min(dy);
end

% f_compute_tke_budget
function [T P D S C] = f_compute_tke_budget(dvel_ds, dvelmean_dy, dvelfluc_ds, dpmean_dy, rho, vel_fluc, pres_fluc, mth, y_norm_mth, timestep)

    load variables/global_parameters.mat;
    load variables/grid.mat;

    % Compute strain-rate tensor
    [SS S11 S12 S13 S22 S23 S33] = f_compute_strain_rate_tensor(dvel_ds);
    S21 = S12; S31 = S13; S32 = S23;

    % Viscous stress tensor
    divU = squeeze(dvel_ds(1,1,:,:,:) + dvel_ds(2,2,:,:,:) + dvel_ds(3,3,:,:,:));
    sig = zeros(3,3,mp,np,pp);
    sig(1,1,:,:,:) = 2/Re*(S11 - divU/3);
    sig(1,2,:,:,:) = 2/Re*(S12);
    sig(1,3,:,:,:) = 2/Re*(S13);
    sig(2,1,:,:,:) = 2/Re*(S21);
    sig(2,2,:,:,:) = 2/Re*(S22 - divU/3);
    sig(2,3,:,:,:) = 2/Re*(S23);
    sig(3,1,:,:,:) = 2/Re*(S31);
    sig(3,2,:,:,:) = 2/Re*(S32);
    sig(3,3,:,:,:) = 2/Re*(S33 - divU/3);

    sigmean = zeros(3,3,np);
    dsigmean_dy = zeros(3,3,np);
    for i = 1:3
        for j = 1:3
            sigmean(i,j,:) = squeeze(mean(sig(i,j,:,:,:),[3 5]));
            dsigmean_dy(i,j,:) = f_compute_derivative_1d(sigmean(i,j,:),y_cc);
        end
    end
    sigfluc = sig - permute(repmat(sigmean, [1 1 1 mp pp]), [1 2 4 3 5]);
    clear sig SS S11 S12 S13 S21 S22 S23 S31 S32 S33
    
    % Transport (T)
    T1 = -0.5*mean(...
         rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(1,:,:,:).*vel_fluc(2,:,:,:)) ...
       + rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:)) ...
       + rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(3,:,:,:).*vel_fluc(2,:,:,:)),[1 3]);
    T2 = -mean(pres_fluc.*squeeze(vel_fluc(2,:,:,:)),[1 3]);
    T3 = mean(...
         squeeze(sigfluc(1,2,:,:,:)).*squeeze(vel_fluc(1,:,:,:)) ...
       + squeeze(sigfluc(2,2,:,:,:)).*squeeze(vel_fluc(2,:,:,:)) ...
       + squeeze(sigfluc(3,2,:,:,:)).*squeeze(vel_fluc(3,:,:,:)),[1 3]);
    T1 = T1 / (8/mth); % normalization
    T2 = T2 / (8/mth); % normalization
    T3 = T3 / (8/mth); % normalization
    T0 = T1 + T2 + T3;
    T = f_compute_derivative_1d(T0,y_cc); 

    % Turbulent production (P)
    P1 = -mean(rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*squeeze(dvelmean_dy(1,:));
    P2 = -mean(rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*squeeze(dvelmean_dy(2,:));
    P3 = -mean(rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*squeeze(dvelmean_dy(3,:));
    P1 = P1 / (8/mth); % normalization
    P2 = P2 / (8/mth); % normalization
    P3 = P3 / (8/mth); % normalization
    P = P1 + P2 + P3;

    % Dissipation (D)
    D = -mean(squeeze(sigfluc(1,1,:,:,:).*dvelfluc_ds(1,1,:,:,:) ...
                    + sigfluc(2,1,:,:,:).*dvelfluc_ds(2,1,:,:,:) ... 
                    + sigfluc(3,1,:,:,:).*dvelfluc_ds(3,1,:,:,:) ...
                    + sigfluc(1,2,:,:,:).*dvelfluc_ds(1,2,:,:,:) ...
                    + sigfluc(2,2,:,:,:).*dvelfluc_ds(2,2,:,:,:) ...
                    + sigfluc(3,2,:,:,:).*dvelfluc_ds(3,2,:,:,:) ...
                    + sigfluc(1,3,:,:,:).*dvelfluc_ds(1,3,:,:,:) ...
                    + sigfluc(2,3,:,:,:).*dvelfluc_ds(2,3,:,:,:)  ...
                    + sigfluc(3,3,:,:,:).*dvelfluc_ds(3,3,:,:,:)), [1 3]);
    D = D / (8/mth); % normalization

    % Pressure-strain (S)
    S = mean(pres_fluc.*squeeze(dvelfluc_ds(1,1,:,:,:)),[1 3]) ...
      + mean(pres_fluc.*squeeze(dvelfluc_ds(2,2,:,:,:)),[1 3]) ...
      + mean(pres_fluc.*squeeze(dvelfluc_ds(3,3,:,:,:)),[1 3]);
    S = S / (8/mth); % normalization

    % Mass flux coupling (C)
    C = mean(squeeze(vel_fluc(1,:,:,:)),[1 3]).*(squeeze(dsigmean_dy(1,2,:))) ...
      + mean(squeeze(vel_fluc(2,:,:,:)),[1 3]).*(squeeze(dsigmean_dy(2,2,:)) - dpmean_dy) ...
      + mean(squeeze(vel_fluc(3,:,:,:)),[1 3]).*(squeeze(dsigmean_dy(3,2,:)));
    C = C / (8/mth); % normalization

    save(post_stat_dir+"/tke_budget_data/tke_budget_"+string(timestep)+".mat","y_norm_mth","T1","T2","T3","P1","P2","P3","D","mth");

end

% f_average_tke_budget
function f_average_tke_budget()

    load variables/global_parameters.mat;
    load variables/grid.mat;
    load variables/options.mat;

    ybeg = -5; yend = 5; ny = 101;
    y = linspace(ybeg,yend,ny);

    T1_averaged = zeros(ny,1);
    T2_averaged = zeros(ny,1);
    T3_averaged = zeros(ny,1);
    P1_averaged = zeros(ny,1);
    P2_averaged = zeros(ny,1);
    P3_averaged = zeros(ny,1);
    D_averaged = zeros(ny,1);
    
    for q = start_idx:Nfiles
        load(post_stat_dir+"/tke_budget_data/tke_budget_"+string(timesteps(q))+".mat");
        i_start = 1;
        for j = 1:ny
            for i = i_start:length(y_norm_mth) - 1
                if (y_norm_mth(i) <= y(j) && y_norm_mth(i+1) > y(j))
                    T1_averaged(j) = T1_averaged(j) + ((T1(i+1) - T1(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + T1(i))/(Nfiles - start_idx + 1);
                    T2_averaged(j) = T2_averaged(j) + ((T2(i+1) - T2(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + T2(i))/(Nfiles - start_idx + 1);
                    T3_averaged(j) = T3_averaged(j) + ((T3(i+1) - T3(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + T3(i))/(Nfiles - start_idx + 1);
                    P1_averaged(j) = P1_averaged(j) + ((P1(i+1) - P1(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + P1(i))/(Nfiles - start_idx + 1);
                    P2_averaged(j) = P2_averaged(j) + ((P2(i+1) - P2(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + P2(i))/(Nfiles - start_idx + 1);
                    P3_averaged(j) = P3_averaged(j) + ((P3(i+1) - P3(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + P3(i))/(Nfiles - start_idx + 1);
                    D_averaged(j) = D_averaged(j) + ((D(i+1) - D(i))/(y_norm_mth(i+1) - y_norm_mth(i))*(y(j) - y_norm_mth(i)) + D(i))/(Nfiles - start_idx + 1);
                    i_start = i;
                    break;
                end
            end
        end
    end 

    T0_averaged = T1_averaged + T2_averaged + T3_averaged;
    T_averaged = f_compute_derivative_1d(T0_averaged,y*mth); 
    P_averaged = P1_averaged + P2_averaged + P3_averaged;

    %
    f1 = figure("DefaultAxesFontSize",18);
    plot(y,T_averaged,'-bo','LineWidth',2); hold on; grid on;
    plot(y,P_averaged,'-go','LineWidth',2);
    plot(y,D_averaged,'-ro','LineWidth',2);
    % xlim([-5 5]); %ylim([-0.002 0.002]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    saveas(f1, post_stat_dir + "/tke_budget/tke_budget_self_similar","png");
    close(f1);
end

% f_compute_energy_spectrum
function [k E] = f_compute_energy_spectrum(vel_fluc, mth)

    load variables/global_parameters.mat;

    [k, E] = compute_energy_spectrum(squeeze(vel_fluc(1,:,floor(np/2),:)), squeeze(vel_fluc(2,:,floor(np/2),:)), squeeze(vel_fluc(3,:,floor(np/2),:)), Lx, Lz);
    E = E/(2^2*mth);
    k = k*mth;
end

% f_compute_vorticity
function omega = f_compute_vorticity(dvel_ds)

    load variables/global_parameters.mat;

    omega = zeros(3,mp,np,pp);

    omega(1,:,:,:) = dvel_ds(3,2,:,:,:) - dvel_ds(2,3,:,:,:);
    omega(2,:,:,:) = dvel_ds(1,3,:,:,:) - dvel_ds(3,1,:,:,:);
    omega(3,:,:,:) = dvel_ds(2,1,:,:,:) - dvel_ds(1,2,:,:,:);

end

%% HELPER FUNCTIONS
% compute the wall-normal derivative of a discretized function, fun(y)
function dfunds = f_compute_derivative_1d(fun,s)

    dfunds = zeros(size(fun));    % initialize discrete derivative vector

    for j = 1:length(s) % for each index in s
        % if at bottom boundary, use 1-sided derivative
        dfunds(1) = (fun(2) - fun(1)) / (s(2) - s(1));
        % if at top boundary, use 1-sided derivative
        dfunds(end) = (fun(end) - fun(end-1)) / (s(end) - s(end-1));
        % otherwise, if in the interior of the domain, 2-sided derivative
        for i = 2:length(s)-1
            dfunds(i) = (fun(i+1) - fun(i-1)) / (s(i+1) - s(i-1));
        end
    end
end

% compute the wall-normal derivative of a discretized function, fun(y)
function dfunds = f_compute_derivative_3d(fun,s,dir)

    dfunds = zeros(size(fun));    % initialize discrete derivative vector

    if (dir == 1)
        % Forward difference at start
        dfunds(1,:,:) = (fun(2,:,:) - fun(1,:,:)) / (s(2) - s(1));
        % Backward difference at end
        dfunds(end,:,:) = (fun(end,:,:) - fun(end-1,:,:)) / (s(end) - s(end-1));
        % Central difference for interior
        for i = 2:length(s)-1
            dfunds(i,:,:) = (fun(i+1,:,:) - fun(i-1,:,:)) / (s(i+1) - s(i-1));
        end
    elseif (dir == 2)
        dfunds(:,1,:) = (fun(:,2,:) - fun(:,1,:)) / (s(2) - s(1));
        dfunds(:,end,:) = (fun(:,end,:) - fun(:,end-1,:)) / (s(end) - s(end-1));
        for i = 2:length(s)-1
            dfunds(:,i,:) = (fun(:,i+1,:) - fun(:,i-1,:)) / (s(i+1) - s(i-1));
        end
    elseif (dir == 3)
        dfunds(:,:,1) = (fun(:,:,2) - fun(:,:,1)) / (s(2) - s(1));
        dfunds(:,:,end) = (fun(:,:,end) - fun(:,:,end-1)) / (s(end) - s(end-1));
        for i = 2:length(s)-1
            dfunds(:,:,i) = (fun(:,:,i+1) - fun(:,:,i-1)) / (s(i+1) - s(i-1));
        end
    else
        disp("Wrong dir input");
    end
end

%% PLOTS
% plot_pdf_pressure
function plot_pdf_pressure(pres,timestep)

    load variables/multiscale.mat;
    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);
    
    histogram(reshape(pres,[],1),[-3:0.002:1.5],'EdgeColor','k','LineWidth',1.5,'Normalization','pdf', 'DisplayStyle', 'stairs'); hold on; grid on;
    plot([pv pv],[1e-8 1e4],'r--','LineWidth',1.5);
    xlim([-3 2]); 
    set(gca, 'YScale', 'log');
    xlabel('$p$','interpreter','latex');
    ylabel('$PDF$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    saveas(f1,post_stat_dir+"/pdf_pressure/pdf_tstep_"+string(timestep),"png"); 
    close(f1);
end

% plot_pdf_omega_xy
function plot_pdf_omega_xy(omega_xy,timestep)

    load variables/linestyle.mat;
    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);
    histogram(reshape(omega_xy,[],1),'EdgeColor','k','LineWidth',1.5,'Normalization','pdf', 'DisplayStyle', 'stairs'); hold on; grid on;
    xlim([0 5]); 
    set(gca, 'YScale', 'log');
    xlabel('$\omega_{xy}$','interpreter','latex');
    ylabel('$PDF$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    saveas(f1,post_stat_dir + "/pdf_omega_xy/pdf_tstep_"+string(timestep),"png"); 
    close(f1);
end

% plot_jpdf_pres_omega_xy
function plot_jpdf_pres_omega_xy(pres, omega_xy, timestep)

    load variables/linestyle.mat;
    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);

    x = reshape(pres,[],1);
    y = reshape(omega_xy,[],1);

    [counts, xEdges, yEdges] = histcounts2(x, y, 100);

    % Convert histogram counts to probability density
    binWidthX = xEdges(2) - xEdges(1);
    binWidthY = yEdges(2) - yEdges(1);
    jointPDF = counts / (sum(counts(:)) * binWidthX * binWidthY);

    % Define bin centers
    xCenters = xEdges(1:end-1) + binWidthX/2;
    yCenters = yEdges(1:end-1) + binWidthY/2;

    % Plot joint PDF as a contour plot
    contourf(xCenters, yCenters, log(jointPDF'), 20, 'LineColor', 'none'); hold on;
    plot([pv pv],[0 5],'r--','LineWidth',1.5);
    xlim([-3 2]); xticks([-3:1:2]);
    ylim([0 5]); yticks([0:1:5]);
    colorbar; caxis([-10 6]);
    xlabel('$p$','Interpreter','latex');
    ylabel('$\omega_{xy}$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    saveas(f1,post_stat_dir + "/jpdf_pres_omega_xy/jpdf_tstep_"+string(timestep),"png"); 
    close(f1);
end

% plot_energy_spectrum
function plot_energy_spectrum(k, E, timestep)

    load variables/linestyle.mat;
    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);

    loglog(k, E, 'b-', 'LineWidth', 2); hold on; grid on;
    loglog([10^(-0.5) 10^(0.3)], [10^(-3) 10^(-13/3)],'k-','LineWidth', 2);
    xlabel('Wavenumber $k$','interpreter','latex');
    ylabel('Energy $E(k)$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    saveas(f1,post_stat_dir + "/energy_spectrum/energy_spectrum_tstep_"+string(timestep),"png"); 
    close(f1);
end

% plot_tke_budget
function plot_tke_budget(y_norm_mth, T, P, D, S, C, timestep)

    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);

    plot(y_norm_mth,T,'-bo','LineWidth',2); hold on; grid on;
    plot(y_norm_mth,P,'-go','LineWidth',2);
    plot(y_norm_mth,D,'-ro','LineWidth',2);
    % plot(y_norm_mth,S,'-mo','LineWidth',2);
    % plot(y_norm_mth,C,'-go','LineWidth',2);
    xlim([-5 5]);
    % ylim([-0.01 0.02]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
    saveas(f1, post_stat_dir + "/tke_budget/tke_budget_tstep_"+string(timestep),"png");
    close(f1);
end

% plot_tke_budget
function plot_tke_budget_T_components(y_norm_mth, T0, T1, T2, T3, timestep)

    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);
    plot(y_norm_mth',T0,'-k','LineWidth',4); hold on; grid on;
    plot(y_norm_mth',T1,'-b^','LineWidth',2);
    plot(y_norm_mth',T2,'-ro','LineWidth',2);
    plot(y_norm_mth',T3,'-g>','LineWidth',2);
    xlim([-5 5]);
    % ylim([-0.01 0.02]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
    saveas(f1, post_stat_dir + "/tke_budget/T_components_tstep_"+string(timestep),"png");
    close(f1);
end

% plot_tke_budget
function plot_tke_budget_Tderiv_components(y_norm_mth, T0, T1, T2, T3, timestep)

    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);
    plot(y_norm_mth',T0,'-k','LineWidth',4); hold on; grid on;
    plot(y_norm_mth',T1,'-b^','LineWidth',2);
    plot(y_norm_mth',T2,'-ro','LineWidth',2);
    plot(y_norm_mth',T3,'-g>','LineWidth',2);
    xlim([-5 5]);
    % ylim([-0.01 0.02]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
    saveas(f1, post_stat_dir + "/tke_budget/Tderiv_components_tstep_"+string(timestep),"png");
    close(f1);
end

% plot_tke_budget
function plot_tke_budget_P_components(y_norm_mth, P, P1, P2, P3, timestep)

    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);
    plot(y_norm_mth',P,'-k','LineWidth',4); hold on; grid on;
    plot(y_norm_mth',P1,'-b^','LineWidth',2);
    plot(y_norm_mth',P2,'-ro','LineWidth',2);
    plot(y_norm_mth',P3,'-g>','LineWidth',2);
    xlim([-5 5]);
    % ylim([-0.01 0.02]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
    saveas(f1, post_stat_dir + "/tke_budget/P_components_tstep_"+string(timestep),"png");
    close(f1);
end

% plot_min_pressure
function plot_min_pressure(time, pres_min)

    load variables/linestyle.mat;
    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize", 18);     

    plot(time,pres_min,'-ko','LineWidth',2); hold on; grid on;
    plot([0 max(time)], [pv pv], 'r--','LineWidth',2);
    axis([0 max(time) -4 1]);
    yticks([-4:1:1]);
    xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
    ylabel('$p_{\mbox{min}}$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
    saveas(f1, post_stat_dir + "/min_pressure.png");
    close(f1);
end

% plot_mom_thickness
function plot_mom_thickness(time, mth)

    load variables/linestyle.mat;
    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize", 18); 
    
    plot(time,mth,'-ko','LineWidth',2); hold on; grid on;
    axis([min(time) max(time) 0 16]);
    xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
    ylabel('$\delta / \delta_\omega^0$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    saveas(f1, post_stat_dir + "/momentum_thickness.png"); 
    close(f1);
end

% plot_growth_rate
function plot_growth_rate(time, mth)

    load variables/global_parameters.mat;

    % Growth rate
    dmth = (mth(2:end) - mth(1:end-1)) ./ (time(2:end) - time(1:end-1));

    % Integrated growth rate
    f_growth_rate = figure("DefaultAxesFontSize", 18); 

    plot(time(1:end-1),dmth,'-ko','LineWidth',2); hold on; grid on;
    plot([0 max(time)],[0.012 0.012], 'r--','LineWidth',1);
    plot([0 max(time)],[0.0135 0.0135], 'b--','LineWidth',1);
    plot([0 max(time)],[0.014 0.014], 'g--','LineWidth',1);
    plot([0 max(time)],[0.017 0.017], 'm--.','LineWidth',1);
    axis([0 max(time) 0 0.1]);
    xlabel('$t \Delta U / \delta_\theta^0$','interpreter','latex');
    ylabel('$\dot{\delta}_{\theta} / \Delta U$','interpreter','latex');
    legend("$\mbox{Present}$",...
            "$0.012^{[1]}$",...     % [1] Baltzer & Livescu (2020, JFM) 
            "$0.0135^{[2]}$",...    % [2] Vaghefi (2014, Thesis) 
            "$0.014^{[3,4]}$",...   % [3] Rogers & Moser (1993, PoF) [4] Blakeley et al. (2023, JFM)
            "$0.017^{[5]}$",...     % [5] Almagro et al. (2017, JFM)
            'Interpreter','latex','location','northeast');
    set(gca,'TickLabelInterpreter','latex');

    saveas(f_growth_rate, post_stat_dir+"/growth_rate.png"); 
    close(f_growth_rate); 
    disp("f_growth_rate saved");
end

% plot_Reynolds_stress
function plot_Reynolds_stress(y_norm_mth, ruu, rvv, rww, ruv, timestep)

    load variables/linestyle.mat;
    load variables/global_parameters.mat;

    ruu = sqrt(ruu);
    rvv = sqrt(rvv);
    rww = sqrt(rww);
    ruv(ruv > 0) = 0;
    ruv = sqrt(-ruv);

    f1 = figure("DefaultAxesFontSize", 18); 

    % sqrt(ruu)/Delta U
    A = readmatrix("data/Bell_Mehta_1990/ruu.dat");
    y_norm_ref1 = A(:,1); ruu_ref1 = A(:,2);
    A = readmatrix("data/Vaghefi_2014/ruu.dat");
    y_norm_ref2 = A(:,1); ruu_ref2 = A(:,2);
    A = readmatrix("data/Wang_et_al_2022/ruu.dat");
    y_norm_ref3 = A(:,1); ruu_ref3 = A(:,2);
    subplot(2,2,1); hold on;
    plot(y_norm_ref1,ruu_ref1,'g+','LineWidth',2);
    plot(y_norm_ref2,ruu_ref2,'b-.','LineWidth',1);
    plot(y_norm_ref3,ruu_ref3,'r--','LineWidth',1);
    plot(y_norm_mth,ruu,'k-','LineWidth',1); grid on;
    axis([-1.5 1.5 0 0.25]);
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.25]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{\left< \overline{\rho} u^{\prime 2} \right>} / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    % sqrt(rvv)/Delta U
    A = readmatrix("data/Bell_Mehta_1990/rvv.dat");
    y_norm_ref1 = A(:,1); rvv_ref1 = A(:,2);
    A = readmatrix("data/Vaghefi_2014/rvv.dat");
    y_norm_ref2 = A(:,1); rvv_ref2 = A(:,2);
    A = readmatrix("data/Wang_et_al_2022/rvv.dat");
    y_norm_ref3 = A(:,1); rvv_ref3 = A(:,2);
    subplot(2,2,2); hold on;
    plot(y_norm_ref1,rvv_ref1,'g+','LineWidth',2);
    plot(y_norm_ref2,rvv_ref2,'b-.','LineWidth',1);
    plot(y_norm_ref3,rvv_ref3,'r--','LineWidth',1);
    plot(y_norm_mth,rvv,'k-','LineWidth',1); grid on;
    axis([-1.5 1.5 0 0.25]); 
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.25]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{\left< \overline{\rho} v^{\prime 2} \right>} / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    % sqrt(rww)/Delta U
    A = readmatrix("data/Bell_Mehta_1990/rww.dat");
    y_norm_ref1 = A(:,1); rww_ref1 = A(:,2);
    A = readmatrix("data/Vaghefi_2014/rww.dat");
    y_norm_ref2 = A(:,1); rww_ref2 = A(:,2);    
    A = readmatrix("data/Wang_et_al_2022/rww.dat");
    y_norm_ref3 = A(:,1); rww_ref3 = A(:,2);
    subplot(2,2,3); hold on;
    plot(y_norm_ref1,rww_ref1,'g+','LineWidth',2);
    plot(y_norm_ref2,rww_ref2,'b-.','LineWidth',1);
    plot(y_norm_ref3,rww_ref3,'r--','LineWidth',1);
    plot(y_norm_mth,rww,'k-','LineWidth',1); grid on;
    axis([-1.5 1.5 0 0.25]); 
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.25]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{\left< \overline{\rho} w^{\prime 2} \right>} / \Delta U$','Interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    % sqrt(-rvu)/Delta U
    A = readmatrix("data/Bell_Mehta_1990/ruv.dat");
    y_norm_ref1 = A(:,1); ruv_ref1 = A(:,2);
    A = readmatrix("data/Vaghefi_2014/ruv.dat");
    y_norm_ref2 = A(:,1); ruv_ref2 = A(:,2);   
    A = readmatrix("data/Wang_et_al_2022/ruv.dat");
    y_norm_ref3 = A(:,1); ruv_ref3 = A(:,2); 
    subplot(2,2,4); hold on;
    plot(y_norm_ref1,ruv_ref1,'g+','LineWidth',2);
    plot(y_norm_ref2,ruv_ref2,'b-.','LineWidth',1);
    plot(y_norm_ref3,ruv_ref3,'r--','LineWidth',1);
    plot(y_norm_mth,ruv,'k-','LineWidth',1); grid on;
    axis([-1.5 1.5 0 0.25]); 
    xticks([-1.5:0.5:1.5]);
    yticks([0:0.05:0.25]);
    xlabel('$y/\delta_\omega$','Interpreter','latex'); 
    ylabel('$\sqrt{-\left< \overline{\rho} u^{\prime} v^{\prime}\right>}  / \Delta U$','Interpreter','latex');
    legend("$\mbox{Present}$",...
            "$\mbox{Bell & Mehta (1990)}$",...
            "$\mbox{Vaghefi (2014)$",...
            "$\mbox{Wang et al, (2022)}$",...
            'Interpreter','latex','location','northeast');
    set(gca,'TickLabelInterpreter','latex');

    saveas(f1, post_stat_dir + "/Reynolds_stress/Reynolds_stress_"+string(timestep),"png");
    close(f1);

end

% plot_integrated_tke_budget
function plot_integrated_tke_budget(time, Pintegrated, Dintegrated)

    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize", 18); 
    
    plot(time,Pintegrated,'-ro','LineWidth',2); hold on; grid on;
    plot(time,Dintegrated,'-bo','LineWidth',2);
    xlim([min(time) max(time)]);
    xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
    % ylabel('$\delta / \delta_\omega^0$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    saveas(f1, post_stat_dir + "/integrated_tke_budget.png"); 
    close(f1);
end


%% PRINT DATA
% print_kolmogorov_length
function print_kolmogorov_length(time, eta_min)

    load variables/global_parameters.mat;

    fileID = fopen(post_stat_dir + "/kolmogorov_length.dat",'w');
    fprintf(fileID,'%12.8f %12.8f %12.8f\r\n',[time; eta_min]);
    fclose(fileID);

end

% print_min_pressure
function print_min_pressure(time, pres_min)

    load variables/global_parameters.mat;

    fileID = fopen(post_stat_dir + "/min_pressure.dat",'w');
    fprintf(fileID,'%12.8f %12.8f\r\n',[time; pres_min]);
    fclose(fileID);
end

% print_mixlayer_thickness
function print_mixlayer_thickness(time, mth, vth)

    load variables/global_parameters.mat;

    fileID = fopen(post_stat_dir + "/mixlayer_thickness.dat",'w');
    fprintf(fileID,'%12.8f %12.8f %12.8f\r\n',[time; mth; vth]);
    fclose(fileID);
end

