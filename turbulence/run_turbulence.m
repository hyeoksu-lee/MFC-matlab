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

    set_global_parameters(l); load variables/global_parameters.mat;

    set_grid(); load variables/grid.mat;

    set_allocate_temporal_data();

    % Parallized loop over timesteps
    if (read_raw_data == "T")
        for q = 1:Nfiles

            % Read data in conservative form
            qc = f_read_data(strcat(mfc_dir(l),"/restart_data/lustre_",int2str(timesteps(q)),".dat"));

            % Convert conservative variables to primitive variables
            qp = f_convert_cons_to_prim(qc);

            % Compute minimum pressure PDF of minimum pressure
            pres_min(q) = min(qp(E_idx,:,:,:),[],"all");

            % Compute PDF of pressure
            plot_pdf_pressure(qp(E_idx,:,:,:), timesteps(q));

            % Compute mean quantities
            [qp_mean qp_fluc] = f_compute_qp_mean_fluc(qp);
            
            % Save mean streamwise velocity
            u_mean(:,q) = qp_mean(momxb,:)';

            % Compute derivatives
            [dvel_ds dvelmean_dy dvelfluc_ds] = f_compute_vel_derivatives(qp(momxb:momxe,:,:,:), qp_mean(momxb:momxe,:), qp_fluc(momxb:momxe,:,:,:));
            dpmean_dy = f_compute_derivative_1d(squeeze(qp_mean(E_idx,:)),y_cc);

            % Compute vorticity
            omega = f_compute_vorticity(dvel_ds);
            omega_xy = sqrt(omega(1,:,:,:).^2 + omega(2,:,:,:).^2);
            plot_pdf_omega_xy(omega_xy,timesteps(q));
            plot_jpdf_pres_omega_xy(qp(E_idx,:,:,:),omega_xy,timesteps(q));

            % Compute mixing layer thickness
            [vth(q) mth(q) y_norm(:,q)] = f_compute_mixing_layer_thickness(qp_mean,dvelmean_dy);

            % Compute Reynolds stress
            [ruu(:,q) rvv(:,q) rww(:,q) ruv(:,q)] = f_compute_Reynolds_stress(qp(1,:,:,:), qp_fluc(momxb:momxe,:,:,:));
            
            % Compute Kolmogorov length scale / dy_min
            eta_min(q) = f_compute_kolmogorov_scale(dvelfluc_ds);

            % Compute TKE budget
            [T P D S C] = f_compute_tke_budget(dvel_ds, dvelmean_dy, dvelfluc_ds, dpmean_dy, squeeze(qp(1,:,:,:)), qp_fluc(momxb:momxe,:,:,:), squeeze(qp(E_idx,:,:,:)), mth(q), y_norm(:,q), timesteps(q));
            plot_tke_budget(y_norm(:,q), T, P, D, S, C, timesteps(q));

            % Compute 1D energy spectrum
            [k,E] = f_compute_energy_spectrum(qp_fluc, mth(q));
            plot_energy_spectrum(k, E, timesteps(q));

        end

        % SAVE DATA AS .MAT FILES
        save(post_stat_dir+"/mean_streamwise_vel.mat","time","y_norm","u_mean");
        save(post_stat_dir+"/Reynolds_stress.mat","time","y_norm","ruu","rvv","rww","ruv");
        save(post_stat_dir+"/mixlayer_thickness.mat","time","mth","vth");
        save(post_stat_dir+"/tke_budget.mat","y_norm","T","P","D","S","C");

        % SAVE DATA AS .DAT FILES
        print_kolmogorov_length(time, eta_min);
        print_min_pressure(time, pres_min);
        print_mixlayer_thickness(time, mth, vth);

        % PLOT
        plot_min_pressure(time, pres_min);
        plot_mom_thickness(time, mth);
    end

    if (avg_over_self_sim == "T")
        f_average_tke_budget();
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

    % Compute fluctuating quantities
    qp_fluc = zeros(sys_size,mp,np,pp);

    % density fluctuation
    qp_fluc(1,:,:,:) = qp(1,:,:,:) - permute(repmat(qp_mean(1,:), [1 1 mp pp]), [1 3 2 4]);

    % velocity fluctuation
    qp_fluc(momxb:momxe,:,:,:) = qp(momxb:momxe,:,:,:) - permute(repmat(qp_mean(momxb:momxe,:), [1 1 mp pp]), [1 3 2 4]);

    % pressure fluctuatino
    qp_fluc(E_idx,:,:,:) = qp(E_idx,:,:,:) - permute(repmat(qp_mean(E_idx,:,:,:), [1 1 mp pp]), [1 3 2 4]);

end

% f_compute_vel_derivatives
function [dvel_ds dvelmean_dy dvelfluc_ds] = f_compute_vel_derivatives(vel, vel_mean, vel_fluc)

    load variables/global_parameters.mat;
    load variables/grid.mat;
    
    dvel_ds = zeros(3,3,mp,np,pp);
    dvelmean_dy = zeros(3,np);
    dvelfluc_ds = zeros(3,3,mp,np,pp);

    % Compute velocity derivatives
    dvel_ds(1,1,:,:,:) = f_compute_derivative_3d(squeeze(vel(1,:,:,:)),x_cc,1);
    dvel_ds(2,1,:,:,:) = f_compute_derivative_3d(squeeze(vel(2,:,:,:)),x_cc,1);
    dvel_ds(3,1,:,:,:) = f_compute_derivative_3d(squeeze(vel(3,:,:,:)),x_cc,1);

    dvel_ds(1,2,:,:,:) = f_compute_derivative_3d(squeeze(vel(1,:,:,:)),y_cc,2);
    dvel_ds(2,2,:,:,:) = f_compute_derivative_3d(squeeze(vel(2,:,:,:)),y_cc,2);
    dvel_ds(3,2,:,:,:) = f_compute_derivative_3d(squeeze(vel(3,:,:,:)),y_cc,2);

    dvel_ds(1,3,:,:,:) = f_compute_derivative_3d(squeeze(vel(1,:,:,:)),z_cc,3);
    dvel_ds(2,3,:,:,:) = f_compute_derivative_3d(squeeze(vel(2,:,:,:)),z_cc,3);
    dvel_ds(3,3,:,:,:) = f_compute_derivative_3d(squeeze(vel(3,:,:,:)),z_cc,3);

    % favre-averaged velocity derivatives
    dvelmean_dy(1,:) = f_compute_derivative_1d(vel_mean(1,:),y_cc);
    dvelmean_dy(2,:) = f_compute_derivative_1d(vel_mean(2,:),y_cc);
    dvelmean_dy(3,:) = f_compute_derivative_1d(vel_mean(3,:),y_cc);

    % fluctuating velocity derivatives
    dvelfluc_ds(1,1,:,:,:) = f_compute_derivative_3d(squeeze(vel_fluc(1,:,:,:)),x_cc,1);
    dvelfluc_ds(2,1,:,:,:) = f_compute_derivative_3d(squeeze(vel_fluc(2,:,:,:)),x_cc,1);
    dvelfluc_ds(3,1,:,:,:) = f_compute_derivative_3d(squeeze(vel_fluc(3,:,:,:)),x_cc,1);

    dvelfluc_ds(1,2,:,:,:) = f_compute_derivative_3d(squeeze(vel_fluc(1,:,:,:)),y_cc,2);
    dvelfluc_ds(2,2,:,:,:) = f_compute_derivative_3d(squeeze(vel_fluc(2,:,:,:)),y_cc,2);
    dvelfluc_ds(3,2,:,:,:) = f_compute_derivative_3d(squeeze(vel_fluc(3,:,:,:)),y_cc,2);

    dvelfluc_ds(1,3,:,:,:) = f_compute_derivative_3d(squeeze(vel_fluc(1,:,:,:)),z_cc,3);
    dvelfluc_ds(2,3,:,:,:) = f_compute_derivative_3d(squeeze(vel_fluc(2,:,:,:)),z_cc,3);
    dvelfluc_ds(3,3,:,:,:) = f_compute_derivative_3d(squeeze(vel_fluc(3,:,:,:)),z_cc,3);
end

% f_compute_mixing_layer_thickness
function [vth mth y_norm] = f_compute_mixing_layer_thickness(qp_mean, dvelmean_dy)

    load variables/grid.mat;
    load variables/index.mat;

    % Compute vorticity thickness
    vth = 2 / max(abs(dvelmean_dy(1,:)),[],"all");
    
    % Compute momentum thickness
    f = qp_mean(1,:) .* (1 - qp_mean(momxb, :)) .* (1 + qp_mean(momxb, :)) / 4;
    mth = trapz(y_cc,f);

    % Compute mth-normalized y value
    y_norm = y_cc'/mth;

end

% f_compute_Reynolds_stress
function [ruu rvv rww ruv] = f_compute_Reynolds_stress(rho, vel_fluc)
    % Compute Reynolds stress
    ruu = mean(squeeze(rho) .* squeeze(vel_fluc(1, :, :, :)).^2, [1 3]);
    rvv = mean(squeeze(rho) .* squeeze(vel_fluc(2, :, :, :)).^2, [1 3]);
    rww = mean(squeeze(rho) .* squeeze(vel_fluc(3, :, :, :)).^2, [1 3]);
    ruv = mean(squeeze(rho) .* squeeze(vel_fluc(1, :, :, :)) .* squeeze(vel_fluc(2, :, :, :)), [1 3]);
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
function [T P D S C] = f_compute_tke_budget(dvel_ds, dvelmean_dy, dvelfluc_ds, dpmean_dy, rho, vel_fluc, pres_fluc, mth, y_norm, timestep)

    load variables/global_parameters.mat;
    load variables/grid.mat;
    load variables/options.mat;

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
    T0 = T1 + T2 + T3;
    plot_tke_budget_T_components(y_norm, T0, T1, T2, T3, timestep);
    T = f_compute_derivative_1d(T0,y_cc); 
    T1deriv = f_compute_derivative_1d(T1,y_cc); 
    T2deriv = f_compute_derivative_1d(T2,y_cc); 
    T3deriv = f_compute_derivative_1d(T3,y_cc); 
    plot_tke_budget_Tderiv_components(y_norm, T, T1deriv, T2deriv, T3deriv, timestep);
    % clear T0 T1 T2 T3

    % Turbulent production (P)
    P1 = -mean(rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*squeeze(dvelmean_dy(1,:));
    P2 = -mean(rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*squeeze(dvelmean_dy(2,:));
    P3 = -mean(rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*squeeze(dvelmean_dy(3,:));
    P = P1 + P2 + P3;
    plot_tke_budget_P_components(y_norm, P, P1, P2, P3, timestep);
    % clear P1 P2 P3

    % Dissipation (D)
    D = -mean(squeeze(sigfluc(1,1,:,:,:).*dvelfluc_ds(1,1,:,:,:)) ...
            + squeeze(sigfluc(2,1,:,:,:).*dvelfluc_ds(2,1,:,:,:)) ... 
            + squeeze(sigfluc(3,1,:,:,:).*dvelfluc_ds(3,1,:,:,:)) ...
            + squeeze(sigfluc(1,2,:,:,:).*dvelfluc_ds(1,2,:,:,:)) ...
            + squeeze(sigfluc(2,2,:,:,:).*dvelfluc_ds(2,2,:,:,:)) ...
            + squeeze(sigfluc(3,2,:,:,:).*dvelfluc_ds(3,2,:,:,:)) ...
            + squeeze(sigfluc(1,3,:,:,:).*dvelfluc_ds(1,3,:,:,:)) ...
            + squeeze(sigfluc(2,3,:,:,:).*dvelfluc_ds(2,3,:,:,:))  ...
            + squeeze(sigfluc(3,3,:,:,:).*dvelfluc_ds(3,3,:,:,:)), [1 3]);

    % Pressure-strain (S)
    S = mean(pres_fluc.*squeeze(dvelfluc_ds(1,1,:,:,:)),[1 3]) ...
      + mean(pres_fluc.*squeeze(dvelfluc_ds(2,2,:,:,:)),[1 3]) ...
      + mean(pres_fluc.*squeeze(dvelfluc_ds(3,3,:,:,:)),[1 3]);

    % Mass flux coupling (C)
    C = mean(squeeze(vel_fluc(1,:,:,:)),[1 3]).*(squeeze(dsigmean_dy(1,2,:,:,:))) ...
      + mean(squeeze(vel_fluc(2,:,:,:)),[1 3]).*(squeeze(dsigmean_dy(2,2,:,:,:)) - dpmean_dy) ...
      + mean(squeeze(vel_fluc(3,:,:,:)),[1 3]).*(squeeze(dsigmean_dy(3,2,:,:,:)));

    % Normalization
    T = T / (8/mth); P = P / (8/mth); D = D / (8/mth); S = S / (8/mth); C = C / (8/mth);

    if (avg_over_self_sim == "T")
        save(post_stat_dir+"/tke_budget_data/tke_budget_"+string(timestep)+".mat","y_norm","T1","T2","T3","P1","P2","P3","D","mth");
    end
end

% f_average_tke_budget
function f_average_tke_budget()

    load variables/global_parameters.mat;
    load variables/grid.mat;

    ybeg = -5; yend = 5; ny = 101;
    y = linspace(ybeg,yend,ny);

    T1_averaged = zeros(ny,1);
    T2_averaged = zeros(ny,1);
    T3_averaged = zeros(ny,1);
    P1_averaged = zeros(ny,1);
    P2_averaged = zeros(ny,1);
    P3_averaged = zeros(ny,1);
    D_averaged = zeros(ny,1);

    for q = 1:Nfiles
        load(post_stat_dir+"/tke_budget_data/tke_budget_"+string(timesteps(q))+".mat");
        for j = 1:ny
            for i = 1:length(y_norm) - 1
                if (y_norm(i) <= y(j) && y_norm(i+1) > y(j))
                    T1_averaged(j) = T1_averaged(j) + ((T1(i+1) - T1(i))/(y_norm(i+1) - y_norm(i))*(y(j) - y_norm(i)) + T1(i))/Nfiles/(8/mth);
                    T2_averaged(j) = T2_averaged(j) + ((T2(i+1) - T2(i))/(y_norm(i+1) - y_norm(i))*(y(j) - y_norm(i)) + T2(i))/Nfiles/(8/mth);
                    T3_averaged(j) = T3_averaged(j) + ((T3(i+1) - T3(i))/(y_norm(i+1) - y_norm(i))*(y(j) - y_norm(i)) + T3(i))/Nfiles/(8/mth);
                    P1_averaged(j) = P1_averaged(j) + ((P1(i+1) - P1(i))/(y_norm(i+1) - y_norm(i))*(y(j) - y_norm(i)) + P1(i))/Nfiles/(8/mth);
                    P2_averaged(j) = P2_averaged(j) + ((P2(i+1) - P2(i))/(y_norm(i+1) - y_norm(i))*(y(j) - y_norm(i)) + P2(i))/Nfiles/(8/mth);
                    P3_averaged(j) = P3_averaged(j) + ((P3(i+1) - P3(i))/(y_norm(i+1) - y_norm(i))*(y(j) - y_norm(i)) + P3(i))/Nfiles/(8/mth);
                    D_averaged(j) = D_averaged(j) + ((D(i+1) - D(i))/(y_norm(i+1) - y_norm(i))*(y(j) - y_norm(i)) + D(i))/Nfiles;
                    break;
                end
            end
        end
    end

    T0 = T1_averaged + T2_averaged + T3_averaged;
    T = f_compute_derivative_1d(T0,y); 
    P = P1_averaged + P2_averaged + P3_averaged;

    %
    f1 = figure("DefaultAxesFontSize",18);
    plot(y,T0,'-k','LineWidth',4); hold on; grid on;
    plot(y,T1_averaged,'-b^','LineWidth',2);
    plot(y,T2_averaged,'-ro','LineWidth',2);
    plot(y,T3_averaged,'-g>','LineWidth',2);
    xlim([-5 5]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    saveas(f1, post_stat_dir + "/tke_budget/T_components_self_similar","png");
    close(f1);

    %
    f1 = figure("DefaultAxesFontSize",18);
    plot(y,T,'-bo','LineWidth',2); hold on; grid on;
    plot(y,P,'-go','LineWidth',2);
    plot(y,D_averaged,'-ro','LineWidth',2);
    xlim([-5 5]);
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
        if j == 1   % if at bottom boundary, use 1-sided derivative
            dfunds(j) = (fun(j+1) - fun(j))/(s(j+1) - s(j));
        elseif j == length(s)% if at top boundary, use 1-sided derivative
            dfunds(j) = (fun(j) - fun(j-1))/(s(j) - s(j-1));
        else    % otherwise, if in the interior of the domain, 2-sided derivative
            dfunds(j) = (fun(j+1) - fun(j-1))/((s(j+1) - s(j-1)));
        end
    end
end

% compute the wall-normal derivative of a discretized function, fun(y)
function dfunds = f_compute_derivative_3d(fun,s,dir)

    dfunds = zeros(size(fun));    % initialize discrete derivative vector

    if (dir == 1)
        for j = 1:length(s) % for each index in s
            if j == 1   % if at bottom boundary, use 1-sided derivative
                dfunds(j,:,:) = (fun(j+1,:,:) - fun(j,:,:))/(s(j+1) - s(j));
            elseif j == length(s)% if at top boundary, use 1-sided derivative
                dfunds(j,:,:) = (fun(j,:,:) - fun(j-1,:,:))/(s(j) - s(j-1));
            else    % otherwise, if in the interior of the domain, 2-sided derivative
                dfunds(j,:,:) = (fun(j+1,:,:) - fun(j-1,:,:))/((s(j+1) - s(j-1)));
            end
        end
    elseif (dir == 2)
        for j = 1:length(s) % for each index in s
            if j == 1   % if at bottom boundary, use 1-sided derivative
                dfunds(:,j,:) = (fun(:,j+1,:) - fun(:,j,:))/(s(j+1) - s(j));
            elseif j == length(s)% if at top boundary, use 1-sided derivative
                dfunds(:,j,:) = (fun(:,j,:) - fun(:,j-1,:))/(s(j) - s(j-1));
            else    % otherwise, if in the interior of the domain, 2-sided derivative
                dfunds(:,j,:) = (fun(:,j+1,:) - fun(:,j-1,:))/((s(j+1) - s(j-1)));
            end
        end
    elseif (dir == 3)
        for j = 1:length(s) % for each index in s
            if j == 1   % if at bottom boundary, use 1-sided derivative
                dfunds(:,:,j) = (fun(:,:,j+1) - fun(:,:,j))/(s(j+1) - s(j));
            elseif j == length(s)% if at top boundary, use 1-sided derivative
                dfunds(:,:,j) = (fun(:,:,j) - fun(:,:,j-1))/(s(j) - s(j-1));
            else    % otherwise, if in the interior of the domain, 2-sided derivative
                dfunds(:,:,j) = (fun(:,:,j+1) - fun(:,:,j-1))/((s(j+1) - s(j-1)));
            end
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
    plot([pv pv],[1e-6 1e4],'r--','LineWidth',1.5);
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
function plot_tke_budget(y_norm, T, P, D, S, C, timestep)

    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);

    plot(y_norm,T,'-bo','LineWidth',2); hold on; grid on;
    plot(y_norm,P,'-go','LineWidth',2);
    plot(y_norm,D,'-ro','LineWidth',2);
    % plot(y_norm,S,'-mo','LineWidth',2);
    % plot(y_norm,C,'-go','LineWidth',2);
    xlim([-5 5]);
    % ylim([-0.01 0.02]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
    saveas(f1, post_stat_dir + "/tke_budget/tke_budget_tstep_"+string(timestep),"png");
    close(f1);
end

% plot_tke_budget
function plot_tke_budget_T_components(y_norm, T0, T1, T2, T3, timestep)

    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);
    plot(y_norm',T0,'-k','LineWidth',4); hold on; grid on;
    plot(y_norm',T1,'-b^','LineWidth',2);
    plot(y_norm',T2,'-ro','LineWidth',2);
    plot(y_norm',T3,'-g>','LineWidth',2);
    xlim([-5 5]);
    % ylim([-0.01 0.02]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
    saveas(f1, post_stat_dir + "/tke_budget/T_components_tstep_"+string(timestep),"png");
    close(f1);
end

% plot_tke_budget
function plot_tke_budget_Tderiv_components(y_norm, T0, T1, T2, T3, timestep)

    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);
    plot(y_norm',T0,'-k','LineWidth',4); hold on; grid on;
    plot(y_norm',T1,'-b^','LineWidth',2);
    plot(y_norm',T2,'-ro','LineWidth',2);
    plot(y_norm',T3,'-g>','LineWidth',2);
    xlim([-5 5]);
    % ylim([-0.01 0.02]);
    xlabel('$y/\delta_\theta$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');
    
    saveas(f1, post_stat_dir + "/tke_budget/Tderiv_components_tstep_"+string(timestep),"png");
    close(f1);
end

% plot_tke_budget
function plot_tke_budget_P_components(y_norm, P, P1, P2, P3, timestep)

    load variables/global_parameters.mat;

    f1 = figure("DefaultAxesFontSize",18);
    plot(y_norm',P,'-k','LineWidth',4); hold on; grid on;
    plot(y_norm',P1,'-b^','LineWidth',2);
    plot(y_norm',P2,'-ro','LineWidth',2);
    plot(y_norm',P3,'-g>','LineWidth',2);
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
    plot([0 400], [pv pv], 'r--','LineWidth',2);
    axis([0 400 -4 1]);
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
    axis([0 400 0 14]);
    xlabel('$t U_1 / \delta_\omega^0$','interpreter','latex');
    ylabel('$\delta / \delta_\omega^0$','interpreter','latex');
    set(gca,'TickLabelInterpreter','latex');

    saveas(f1, post_stat_dir + "/momentum_thickness.png"); 
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


