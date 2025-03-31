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
        dpmean_dy = compute_derivative_1d(squeeze(qp_mean(E_idx,:)),y_cc);

        % Compute vorticity
        omega = f_compute_vorticity(dvel_ds);
        omega_xy = sqrt(omega(1,:,:,:).^2 + omega(2,:,:,:).^2);
        plot_pdf_omega_xy(omega_xy,timesteps(q));
        plot_jpdf_pres_omega_xy(qp(E_idx,:,:,:),omega_xy,timesteps(q));

        % Compute mixing layer thickness
        [vth(q) mth(q) y_norm(:,q)] = f_compute_mixing_layer_thickness(qp_mean,dvelmean_dy);

        % Compute Reynolds stress
        [ruu(:, q) rvv(:,q) rww(:,q) ruv(:,q)] = f_compute_Reynolds_stress(qp(1,:,:,:), qp_fluc(momxb:momxe,:,:,:));
        
        % Compute Kolmogorov length scale / dy_min
        eta_min(q) = f_compute_kolmogorov_scale(dvelfluc_ds);

        if (q == Nfiles(l))
            % Compute TKE budget
            [T P D S C] = f_compute_tke_budget(dvel_ds);
            plot_tke_budget(y_norm(:,q),T,P,D,S,C,timesteps(q));

            % Compute 1D energy spectrum
            [k,E] = f_compute_energy_spectrum(qp_fluc);
            plot_energy_spectrum(k,E,timesteps(q));
        end
    end

    % SAVE DATA AS .MAT FILES
    disp('save data')
    save(post_stat_dir+"/mean_streamwise_vel.mat","time","y_norm","u_mean");
    save(post_stat_dir+"/Reynolds_stress.mat","time","y_norm","ruu","rvv","rww","ruv");
    save(post_stat_dir+"/mix_layer_thickness.mat","time","mth","vth");
    save(post_stat_dir+"/tke_budget.mat","y_norm","T","P","D","S","C");

    % SAVE DATA AS .DAT FILES
    print_kolmogorov_length(time, eta_max, eta_min);
    print_min_pressure(time, pres_min);
    print_mixlayer_thickness(time, mth, vth);

    % PLOT
    plot_min_pressure(time, pres_min);
    plot_mom_thickness(time, mth);
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
    qc = reshape(A, sys_size, mp, np, pp);
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
    dvel_ds(1,1,:,:,:) = compute_derivative_3d(squeeze(vel(1,:,:,:)),x_cc,1);
    dvel_ds(2,1,:,:,:) = compute_derivative_3d(squeeze(vel(2,:,:,:)),x_cc,1);
    dvel_ds(3,1,:,:,:) = compute_derivative_3d(squeeze(vel(3,:,:,:)),x_cc,1);

    dvel_ds(1,2,:,:,:) = compute_derivative_3d(squeeze(vel(1,:,:,:)),y_cc,2);
    dvel_ds(2,2,:,:,:) = compute_derivative_3d(squeeze(vel(2,:,:,:)),y_cc,2);
    dvel_ds(3,2,:,:,:) = compute_derivative_3d(squeeze(vel(3,:,:,:)),y_cc,2);

    dvel_ds(1,3,:,:,:) = compute_derivative_3d(squeeze(vel(1,:,:,:)),z_cc,3);
    dvel_ds(2,3,:,:,:) = compute_derivative_3d(squeeze(vel(2,:,:,:)),z_cc,3);
    dvel_ds(3,3,:,:,:) = compute_derivative_3d(squeeze(vel(3,:,:,:)),z_cc,3);

    % favre-averaged velocity derivatives
    dvelmean_dy(1,:) = compute_derivative_1d(vel_mean(1,:),y_cc);
    dvelmean_dy(2,:) = compute_derivative_1d(vel_mean(2,:),y_cc);
    dvelmean_dy(3,:) = compute_derivative_1d(vel_mean(3,:),y_cc);

    % fluctuating velocity derivatives
    dvelfluc_ds(1,1,:,:,:) = compute_derivative_3d(squeeze(vel_fluc(1,:,:,:)),x_cc,1);
    dvelfluc_ds(2,1,:,:,:) = compute_derivative_3d(squeeze(vel_fluc(2,:,:,:)),x_cc,1);
    dvelfluc_ds(3,1,:,:,:) = compute_derivative_3d(squeeze(vel_fluc(3,:,:,:)),x_cc,1);

    dvelfluc_ds(1,2,:,:,:) = compute_derivative_3d(squeeze(vel_fluc(1,:,:,:)),y_cc,2);
    dvelfluc_ds(2,2,:,:,:) = compute_derivative_3d(squeeze(vel_fluc(2,:,:,:)),y_cc,2);
    dvelfluc_ds(3,2,:,:,:) = compute_derivative_3d(squeeze(vel_fluc(3,:,:,:)),y_cc,2);

    dvelfluc_ds(1,3,:,:,:) = compute_derivative_3d(squeeze(vel_fluc(1,:,:,:)),z_cc,3);
    dvelfluc_ds(2,3,:,:,:) = compute_derivative_3d(squeeze(vel_fluc(2,:,:,:)),z_cc,3);
    dvelfluc_ds(3,3,:,:,:) = compute_derivative_3d(squeeze(vel_fluc(3,:,:,:)),z_cc,3);
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
function [T P D S C] = f_compute_tke_budget(dvel_ds)

    load variables/global_parameters.mat;

    % Compute strain-rate tensor
    [SS S11 S12 S13 S22 S23 S33] = f_compute_strain_rate_tensor(dvel_ds);
    S21 = S12; S31 = S13; S32 = S23;

    % Viscous stress tensor
    divU = squeeze(dvel_ds(1,1,:,:,:) + dvel_ds(2,2,:,:,:) + dvel_ds(3,3,:,:,:));
    v11 = 2/Re*(S11 - divU/3); v12 = 2/Re*(S12); v13 = 2/Re*(S13);
    v21 = 2/Re*(S21); v22 = 2/Re*(S22 - divU/3); v23 = 2/Re*(S23);
    v31 = 2/Re*(S31); v32 = 2/Re*(S32); v33 = 2/Re*(S33 - divU/3);
    v11_mean = mean(v11,[1 3]);
    v12_mean = mean(v12,[1 3]);
    v13_mean = mean(v13,[1 3]);
    v21_mean = mean(v21,[1 3]);
    v22_mean = mean(v22,[1 3]);
    v23_mean = mean(v23,[1 3]);
    v31_mean = mean(v31,[1 3]);
    v32_mean = mean(v32,[1 3]);
    v33_mean = mean(v33,[1 3]);
    
    dv11_meandx = compute_derivative_1d(v11_mean,x_cc);
    dv12_meandy = compute_derivative_1d(v12_mean,y_cc);
    dv13_meandz = compute_derivative_1d(v13_mean,z_cc);
    dv21_meandx = compute_derivative_1d(v21_mean,x_cc);
    dv22_meandy = compute_derivative_1d(v22_mean,y_cc);
    dv23_meandz = compute_derivative_1d(v23_mean,z_cc);
    dv31_meandx = compute_derivative_1d(v31_mean,x_cc);
    dv32_meandy = compute_derivative_1d(v32_mean,y_cc);
    dv33_meandz = compute_derivative_1d(v33_mean,z_cc);

    v11_fluc = v11 - squeeze(permute(repmat(v11_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]));
    v12_fluc = v12 - squeeze(permute(repmat(v12_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]));
    v13_fluc = v13 - squeeze(permute(repmat(v13_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]));
    v21_fluc = v21 - squeeze(permute(repmat(v21_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]));
    v22_fluc = v22 - squeeze(permute(repmat(v22_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]));
    v23_fluc = v23 - squeeze(permute(repmat(v23_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]));
    v31_fluc = v31 - squeeze(permute(repmat(v31_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]));
    v32_fluc = v32 - squeeze(permute(repmat(v32_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]));
    v33_fluc = v33 - squeeze(permute(repmat(v33_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]));
    
    % % Transport (T)
    % T1 = rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(1,:,:,:).*vel_fluc(1,:,:,:));
    % T2 = rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:).*vel_fluc(1,:,:,:));
    % T3 = rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(3,:,:,:).*vel_fluc(1,:,:,:));
    % T4 = mean(pres_fluc.*squeeze(vel_fluc(1,:,:,:)),[1 3]);
    % T5 = -squeeze(v11_fluc.*squeeze(vel_fluc(1,:,:,:)));
    % T6 = -squeeze(v21_fluc.*squeeze(vel_fluc(2,:,:,:)));
    % T7 = -squeeze(v31_fluc.*squeeze(vel_fluc(3,:,:,:)));
    % T0 = -(0.5*mean(T1 + T2 + T3,[1 3]) + T4 + mean(T5 + T6 + T7,[1 3]));
    % Tx = compute_derivative_1d(T0,x);

    % T1 = rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(1,:,:,:).*vel_fluc(2,:,:,:));
    % T2 = rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:));
    % T3 = rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(3,:,:,:).*vel_fluc(2,:,:,:));
    % T4 = mean(pres_fluc.*squeeze(vel_fluc(2,:,:,:)),[1 3]);
    % T5 = -squeeze(v12_fluc.*squeeze(vel_fluc(1,:,:,:)));
    % T6 = -squeeze(v22_fluc.*squeeze(vel_fluc(2,:,:,:)));
    % T7 = -squeeze(v32_fluc.*squeeze(vel_fluc(3,:,:,:)));
    % T0 = -(0.5*mean(T1 + T2 + T3,[1 3]) + T4 + mean(T5 + T6 + T7,[1 3]));
    % Ty = compute_derivative_1d(T0,y);

    % T1 = rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(1,:,:,:).*vel_fluc(3,:,:,:));
    % T2 = rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:).*vel_fluc(3,:,:,:));
    % T3 = rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(3,:,:,:).*vel_fluc(3,:,:,:));
    % T4 = mean(pres_fluc.*squeeze(vel_fluc(3,:,:,:)),[1 3]);
    % T5 = -squeeze(v13_fluc.*squeeze(vel_fluc(1,:,:,:)));
    % T6 = -squeeze(v23_fluc.*squeeze(vel_fluc(2,:,:,:)));
    % T7 = -squeeze(v33_fluc.*squeeze(vel_fluc(3,:,:,:)));
    % T0 = -(0.5*mean(T1 + T2 + T3,[1 3]) + T4 + mean(T5 + T6 + T7, [1 3]));
    % Tz = compute_derivative_1d(T0,z);

    % T = Tx + Ty + Tz;
    % clear T0 T1 T2 T3 T4 T5 T6 T7 Tx Ty Tz

    % % Turbulent production (P)
    % P1 = mean(rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(1,:,:,:)),[1 3]).*dUmeandx;
    % P2 = mean(rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(1,:,:,:)),[1 3]).*dVmeandx;
    % P3 = mean(rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(1,:,:,:)),[1 3]).*dWmeandx;
    % P4 = mean(rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*dUmeandy;
    % P5 = mean(rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*dVmeandy;
    % P6 = mean(rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(2,:,:,:)),[1 3]).*dWmeandy;
    % P7 = mean(rho.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(3,:,:,:)),[1 3]).*dUmeandz;
    % P8 = mean(rho.*squeeze(vel_fluc(2,:,:,:).*vel_fluc(3,:,:,:)),[1 3]).*dVmeandz;
    % P9 = mean(rho.*squeeze(vel_fluc(3,:,:,:).*vel_fluc(3,:,:,:)),[1 3]).*dWmeandz;
    % P = - (P1 + P2 + P3 + P4 + P5 + P6 + P7 + P8 + P9);
    % clear P1 P2 P3 P4 P5 P6 P7 P8 P9

    % % Dissipation (D)
    % D = -mean(squeeze(v11_fluc).*dudx + squeeze(v21_fluc).*dvdx + squeeze(v31_fluc).*dwdx ...
    %         + squeeze(v12_fluc).*dudy + squeeze(v22_fluc).*dvdy + squeeze(v32_fluc).*dwdy ...
    %         + squeeze(v13_fluc).*dudz + squeeze(v23_fluc).*dvdz + squeeze(v33_fluc).*dwdz, [1 3]);

    % % Pressure-strain (S)
    % S = mean(pres_fluc.*dudx,[1 3]) ...
    %   + mean(pres_fluc.*dvdy,[1 3]) ...
    %   + mean(pres_fluc.*dwdz,[1 3]);

    % % Mass flux coupling (C)
    % C = mean(squeeze(vel_fluc(1,:,:,:)),[1 3]).*(dv11_meandx + dv12_meandy + dv13_meandz - dpmeandx) ...
    %   + mean(squeeze(vel_fluc(2,:,:,:)),[1 3]).*(dv21_meandx + dv22_meandy + dv23_meandz - dpmeandy) ...
    %   + mean(squeeze(vel_fluc(3,:,:,:)),[1 3]).*(dv31_meandx + dv32_meandy + dv33_meandz - dpmeandz);

    % compute production
    dUmeandy = getDerivative(Umean,y_cc); % mean shear
    uvmean = mean(u.*v,[1,3]).'; % mean Reynolds stresses
    P = -uvmean.*dUmeandy; % production term in TKE eqn.

    % compute turbulent transport
    Kpoint = 0.5*(u.^2 + v.^2 + w.^2); % pointwise kinetic energy
    Kpvmean = mean(Kpoint.*v,[1,3]); % turb. wall-normal advection
    T = -getDerivative(Kpvmean,y_cc); % turb. transport in TKE eqn.

    % compute viscous transport
    Kpmean = mean(Kpoint,[1,3]); % wall-normal TKE profile
    viscTrans = nu*getDerivative(getDerivative(Kpmean,y_cc),y_cc); % viscous transport in TKE eqn.
    
    % compute pressure transport
    vpMean = mean(v.*p,[1,3]); % pressure-(v-velocity) corr.
    pressTrans = -getDerivative(vpMean,y_cc); % pressure transport in TKE eqn.
    pressTransPlus = pressTrans/((utau^4)/nu); % normalize in plus units
    
    % compute dissipation
    % call the function provided to the compute strain rate tensor:
    [S11,S12,S13,S22,S23,S33] = compute_Sij(u,v,w,x_cc,y_cc,z_cc);
    Sijsum = S11.^2 + S12.^2 + S13.^2 +... % squared Frobenius norm sum
    S12.^2 + S22.^2 + S23.^2 +...
    S13.^2 + S23.^2 + S33.^2;
    Sijsummean = squeeze(mean(Sijsum,[1,3])).'; % average sum over x and z
    D = -2*nu*Sijsummean; % dissipation in TKE eqn.

end

% f_compute_energy_spectrum
function [k E] = f_compute_energy_spectrum(vel_fluc)
    [k, E] = compute_energy_spectrum(squeeze(vel_fluc(1,:,floor(np(l)/2),:)), squeeze(vel_fluc(2,:,floor(np(l)/2),:)), squeeze(vel_fluc(3,:,floor(np(l)/2),:)), Lx(l), Lz(l));
    E = E/(2^2*mth(q));
    k = k*mth(q);
end

% f_compute_vorticity
function omega = f_compute_vorticity(dvel_ds)

    load variables/global_parameters.mat;

    omega = zeros(3,mp,np,pp);

    omega(1,:,:,:) = dvel_ds(3,2,:,:,:) - dvel_ds(2,3,:,:,:);
    omega(2,:,:,:) = dvel_ds(1,3,:,:,:) - dvel_ds(3,1,:,:,:);
    omega(3,:,:,:) = dvel_ds(2,1,:,:,:) - dvel_ds(1,2,:,:,:);

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

    saveas(f1,post_stat_dir + "/jpdf_pres_vor/jpdf_tstep_"+string(timestep),"png"); 
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
function plot_tke_budget(y_norm,T,P,D,S,C,timestep)

    load variables/global_parameters.mat;

end

% plot_min_pressure
function plot_min_pressure(time, pres_min)

    load variables/linestyle.mat;
    
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
function print_kolmogorov_length(time, eta_max, eta_min)

    load variables/global_parameters.mat;

    fileID = fopen(post_stat_dir + "/kolmogorov_length.dat",'w');
    fprintf(fileID,'%12.8f %12.8f %12.8f\r\n',[time; eta_max; eta_min]);
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


