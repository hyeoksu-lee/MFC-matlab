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
    if (mom_thickness == "T" || mean_streamwise_vel == "T")
        vth = zeros(1,Nfiles(l));       % Vorticity thickness
        mth = zeros(1,Nfiles(l));       % Momentum thickness
    end
    if (Reynolds_stress == "T")
        ruu = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho u u)
        rvv = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho v v)
        rww = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho w w)
        ruv = zeros(np(l),Nfiles(l));   % Reynolds stress: (rho u v)
        eta_min = zeros(1,Nfiles(l));
        eta_max = zeros(1,Nfiles(l));
    end
    if (mean_streamwise_vel == "T")
        u_mean = zeros(np(l),Nfiles(l));% Mean streamwise velocity
        y_norm = zeros(np(l),Nfiles(l));
    end
    if (min_pressure == "T")
        pres_min = zeros(1,Nfiles(l));
    end 

    % Parallized loop over timesteps
    for q = 1:Nfiles(l)

        % Save global parameters
        f_save_global_parameters(l,q);

        % Read data
        qc = f_read_data();


        pres_min(q) = min(pres,[],"all");

        f_pdf_pressure = figure("DefaultAxesFontSize",18); plot_pdf_pressure(f_pdf_pressure, pres);
        saveas(f_pdf_pressure,f_pdf_pressure_dir(l)+string((q-1)*Nt_save(l)),"png"); close(f_pdf_pressure);

        clear A

        if (mom_thickness == "T" || mean_streamwise_vel == "T" || Reynolds_stress == "T" || vortex == "T")
            % Compute mean and fluctuating quantities
            rho_mean = mean(rho,[1 3]); % Mean density
            mom_mean = squeeze(mean(mom,[2 4])); % Mean momentum
            vel_mean = mom_mean./rho_mean; % Favre-averaged velocity
            pres_mean = mean(pres, [1 3]);
        end 

        if (mom_thickness == "T" || mean_streamwise_vel == "T")
            % Compute vorticity thickness
            dUmeandy = compute_derivative_1d(vel_mean(1,:),y);
            dUmeandymax = max(abs(dUmeandy),[],"all");
            vth(q) = 2 / dUmeandymax;
            clear dUmeandy dUmeandymax
            
            % Compute momentum thickness
            f = rho_mean .* (1 - vel_mean(1, :)) .* (1 + vel_mean(1, :)) / 4;
            mth(q) = trapz(y(l,:),f);
        end 

        if (mean_streamwise_vel == "T")
            % Normalized y-axis 
            y_norm(:,q) = y(l,1:np(l))'/mth(q);

            % Mean streamwise velocity
            u_mean(:,q) = vel_mean(1,:)';
        end 

        if (Reynolds_stress == "T")
            % Compute velocity and pressure fluctuation
            vel_fluc = vel - permute(repmat(vel_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]);

            % Compute Reynolds stress
            ruu(:, q) = mean(rho .* squeeze(vel_fluc(1, :, :, :)).^2, [1 3]);
            rvv(:, q) = mean(rho .* squeeze(vel_fluc(2, :, :, :)).^2, [1 3]);
            rww(:, q) = mean(rho .* squeeze(vel_fluc(3, :, :, :)).^2, [1 3]);
            ruv(:, q) = mean(rho .* squeeze(vel_fluc(1, :, :, :)) .* squeeze(vel_fluc(2, :, :, :)), [1 3]);

            % Velocity derivatives
            dUdx = compute_derivative_3d(squeeze(vel(1,:,:,:)),x,1);
            dVdx = compute_derivative_3d(squeeze(vel(2,:,:,:)),x,1);
            dWdx = compute_derivative_3d(squeeze(vel(3,:,:,:)),x,1);

            dUdy = compute_derivative_3d(squeeze(vel(1,:,:,:)),y,2);
            dVdy = compute_derivative_3d(squeeze(vel(2,:,:,:)),y,2);
            dWdy = compute_derivative_3d(squeeze(vel(3,:,:,:)),y,2);

            dUdz = compute_derivative_3d(squeeze(vel(1,:,:,:)),z,3);
            dVdz = compute_derivative_3d(squeeze(vel(2,:,:,:)),z,3);
            dWdz = compute_derivative_3d(squeeze(vel(3,:,:,:)),z,3);

            % favre-averaged velocity derivatives
            dUmeandx = compute_derivative_1d(vel_mean(1,:),x);
            dVmeandx = compute_derivative_1d(vel_mean(2,:),x);
            dWmeandx = compute_derivative_1d(vel_mean(3,:),x);

            dUmeandy = compute_derivative_1d(vel_mean(1,:),y);
            dVmeandy = compute_derivative_1d(vel_mean(2,:),y);
            dWmeandy = compute_derivative_1d(vel_mean(3,:),y);

            dUmeandz = compute_derivative_1d(vel_mean(1,:),z);
            dVmeandz = compute_derivative_1d(vel_mean(2,:),z);
            dWmeandz = compute_derivative_1d(vel_mean(3,:),z);

            % fluctuating velocity derivatives
            dudx = compute_derivative_3d(squeeze(vel_fluc(1,:,:,:)),x,1);
            dvdx = compute_derivative_3d(squeeze(vel_fluc(2,:,:,:)),x,1);
            dwdx = compute_derivative_3d(squeeze(vel_fluc(3,:,:,:)),x,1);

            dudy = compute_derivative_3d(squeeze(vel_fluc(1,:,:,:)),y,2);
            dvdy = compute_derivative_3d(squeeze(vel_fluc(2,:,:,:)),y,2);
            dwdy = compute_derivative_3d(squeeze(vel_fluc(3,:,:,:)),y,2);

            dudz = compute_derivative_3d(squeeze(vel_fluc(1,:,:,:)),z,3);
            dvdz = compute_derivative_3d(squeeze(vel_fluc(2,:,:,:)),z,3);
            dwdz = compute_derivative_3d(squeeze(vel_fluc(3,:,:,:)),z,3);

            % Strain-rate tensor
            S11 = dUdx;
            S12 = 0.5*(dUdy + dVdx);
            S13 = 0.5*(dUdz + dWdx);

            S21 = S12;
            S22 = dVdy;
            S23 = 0.5*(dVdz + dWdy);

            S31 = S13;
            S32 = S23;
            S33 = dWdz;

            SS = S11.^2 + S12.^2 + S13.^2 + S21.^2 + S22.^2 + S23.^2 + S31.^2 + S32.^2 + S33.^2;
            
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
            
            % Kolmogorov length scale / dy_min
            eta_max(q) = (1/(Re0^3 * min(eps_mean)))^(0.25) /dy(l);
            eta_min(q) = (1/(Re0^3 * max(eps_mean)))^(0.25) /dy(l);

            % Pressure statistics
            pres_fluc = pres - squeeze(permute(repmat(pres_mean, [1 1 mp(l) pp(l)]), [1 3 2 4]));
            dpmeandx = compute_derivative_1d(pres_mean,x);
            dpmeandy = compute_derivative_1d(pres_mean,y);
            dpmeandz = compute_derivative_1d(pres_mean,z);

            if (tke_budget == "T" && q == Nfiles(l))
                % Viscous stress tensor
                divU = dUdx + dVdy + dWdz;
                v11 = 2/Re0*(S11 - divU/3); v12 = 2/Re0*(S12); v13 = 2/Re0*(S13);
                v21 = 2/Re0*(S21); v22 = 2/Re0*(S22 - divU/3); v23 = 2/Re0*(S23);
                v31 = 2/Re0*(S31); v32 = 2/Re0*(S32); v33 = 2/Re0*(S33 - divU/3);
                v11_mean = mean(v11,[1 3]);
                v12_mean = mean(v12,[1 3]);
                v13_mean = mean(v13,[1 3]);
                v21_mean = mean(v21,[1 3]);
                v22_mean = mean(v22,[1 3]);
                v23_mean = mean(v23,[1 3]);
                v31_mean = mean(v31,[1 3]);
                v32_mean = mean(v32,[1 3]);
                v33_mean = mean(v33,[1 3]);
                
                dv11_meandx = compute_derivative_1d(v11_mean,x);
                dv12_meandy = compute_derivative_1d(v12_mean,y);
                dv13_meandz = compute_derivative_1d(v13_mean,z);
                dv21_meandx = compute_derivative_1d(v21_mean,x);
                dv22_meandy = compute_derivative_1d(v22_mean,y);
                dv23_meandz = compute_derivative_1d(v23_mean,z);
                dv31_meandx = compute_derivative_1d(v31_mean,x);
                dv32_meandy = compute_derivative_1d(v32_mean,y);
                dv33_meandz = compute_derivative_1d(v33_mean,z);

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
                dUmeandy = getDerivative(Umean,y); % mean shear
                uvmean = mean(u.*v,[1,3]).'; % mean Reynolds stresses
                P = -uvmean.*dUmeandy; % production term in TKE eqn.

                % compute turbulent transport
                Kpoint = 0.5*(u.^2 + v.^2 + w.^2); % pointwise kinetic energy
                Kpvmean = mean(Kpoint.*v,[1,3]); % turb. wall-normal advection
                T = -getDerivative(Kpvmean,y); % turb. transport in TKE eqn.

                % compute viscous transport
                Kpmean = mean(Kpoint,[1,3]); % wall-normal TKE profile
                viscTrans = nu*getDerivative(getDerivative(Kpmean,y),y); % viscous transport in TKE eqn.
                
                % compute pressure transport
                vpMean = mean(v.*p,[1,3]); % pressure-(v-velocity) corr.
                pressTrans = -getDerivative(vpMean,y); % pressure transport in TKE eqn.
                pressTransPlus = pressTrans/((utau^4)/nu); % normalize in plus units
                
                % compute dissipation
                % call the function provided to the compute strain rate tensor:
                [S11,S12,S13,S22,S23,S33] = compute_Sij(u,v,w,x,y,z);
                Sijsum = S11.^2 + S12.^2 + S13.^2 +... % squared Frobenius norm sum
                S12.^2 + S22.^2 + S23.^2 +...
                S13.^2 + S23.^2 + S33.^2;
                Sijsummean = squeeze(mean(Sijsum,[1,3])).'; % average sum over x and z
                D = -2*nu*Sijsummean; % dissipation in TKE eqn.

                save_tke_budget(post_stat_dir(l),y_norm,T,P,D,S,C)
            end

            if (energy_spectrum == "T")
                [k, E] = compute_energy_spectrum(squeeze(vel_fluc(1,:,floor(np(l)/2),:)), squeeze(vel_fluc(2,:,floor(np(l)/2),:)), squeeze(vel_fluc(3,:,floor(np(l)/2),:)), Lx(l), Lz(l));
                E = E/(2^2*mth(q));
                k = k*mth(q);
                f_energy_spectrum = figure("DefaultAxesFontSize",18); plot_energy_spectrum(f_energy_spectrum, k, E);
                saveas(f_energy_spectrum,f_energy_spectrum_dir(l)+string((q-1)*Nt_save(l)),"png"); close(f_energy_spectrum);
            end 

            clear dudx dvdx dwdx dudy dvdy dwdy dudz dvdz dwdz
            clear dUmeandx dVmeandx dWmeandx dUmeandy dVmeandy dWmeandy dUmeandz dVmeandz dWmeandz
            clear v11_fluc v12_fluc v13_fluc v21_fluc v22_fluc v23_fluc v31_fluc v32_fluc v33_fluc

        end

        if (vortex == "T")
            omega_x = dWdy - dVdz;
            omega_y = dUdz - dWdx;
            omega_z = dVdx - dUdy;

            omega_xy = sqrt(omega_x.^2 + omega_y.^2);

            clear dUdx dVdx dWdx dUdy dVdy dWdy dUdz dVdz dWdz

            f_pdf_omega_xy = figure("DefaultAxesFontSize",18); plot_pdf_omega_xy(f_pdf_omega_xy, omega_xy);
            saveas(f_pdf_omega_xy,f_pdf_omega_xy_dir(l)+string((q-1)*Nt_save(l)),"png"); close(f_pdf_omega_xy);

            f_jpdf_pres_vor = figure("DefaultAxesFontSize",18); plot_jpdf_pres_vor(f_jpdf_pres_vor, pres, omega_xy);
            saveas(f_jpdf_pres_vor,f_jpdf_pres_vor_dir(l)+string((q-1)*Nt_save(l)),"png"); close(f_jpdf_pres_vor);
        end
    end

    if (mean_streamwise_vel == "T")
        save_mean_streamwise_vel(post_stat_dir(l), time(l,1:Nfiles(l)), y_norm, u_mean); disp('mean_streamwise_vel data saved');
    end
    if (Reynolds_stress == "T")
        save_Reynolds_stress(post_stat_dir(l), time(l,1:Nfiles(l)), vth, mth, ruu, rvv, rww, ruv); disp('Reynolds_stress data saved');
        print_kolmogorov_length(p_kolmogorov_length_dir(l), time(l,1:Nfiles(l)), eta_max, eta_min); disp('p_kolmogorov_length printed');
    end
    if (min_pressure == "T")
        f_min_pressure = figure("DefaultAxesFontSize", 18); plot_min_pressure(f_min_pressure, time(l,1:Nfiles(l)), pres_min);
        saveas(f_min_pressure, f_min_pressure_dir(l)); close(f_min_pressure); disp('f_min_pressure saved');
        print_min_pressure(p_min_pressure_dir(l), time(l,1:Nfiles(l)), pres_min); disp('p_min_pressure printed');
    end
    if (mom_thickness == "T")
        % Momentum thickness
        f_mom_thickness = figure("DefaultAxesFontSize", 18); plot_mom_thickness(f_mom_thickness, time(l,1:Nfiles(l)), mth);
        saveas(f_mom_thickness, f_mom_thickness_dir(l)); close(f_mom_thickness); disp('f_mom_thickness saved');
        save_mom_thickness(post_stat_dir(l), time(l,1:Nfiles(l)), mth); disp('mom_thickness data saved');
        print_mom_thickness(p_mom_thickness_dir(l), time(l,1:Nfiles(l)), mth); disp('p_momentum_thickness printed');

        % Vorticity thickness
        print_vor_thickness(p_vor_thickness_dir(l), time(l,1:Nfiles(l)), vth); disp('p_vorticity_thickness printed');
    end

end


%% FUNCTIONS
% f_save_global_parameters
function f_save_global_parameters(l, q)
    
    filename = strcat(mfc_dir(l),"/restart_data/lustre_",int2str(timesteps(l,q)),".dat");
    mp = m(l) + 1;
    np = n(l) + 1;
    pp = p(l) + 1;
    r = mp*np*pp;
end

% f_read_data
function qc = f_read_data(filename)

    % Read data
    disp(filename);
    fileID = fopen(filename,'r');
    A = fread(fileID,'double');
    fclose(fileID);

    % Reassign density & velocity components
    r = mp*np*pp;
    rho = reshape(A(1:r), mp, np, pp); % Continuity
    mom = permute(reshape(A(r+1:4*r), mp, np, pp, 3), [4 1 2 3]) * (u_b/u_m); % Momentum
    vel = mom ./ repmat(reshape(rho, 1, mp, np, pp), [3 1 1 1]); % Velocity
    
    E = reshape(A((E_idx - 1)*r + 1:E_idx*r), mp, np, pp) * (u_b/u_m)^2; % Energy
    dyn_p = 0.5 * squeeze(mom(1,:,:,:).*vel(1,:,:,:) + mom(2,:,:,:).*vel(2,:,:,:) + mom(3,:,:,:).*vel(3,:,:,:));
    pres = (1/Gamma) * (E - dyn_p - Pi_inf);

end

% f_
function qp = f_cons_to_prim(qc)