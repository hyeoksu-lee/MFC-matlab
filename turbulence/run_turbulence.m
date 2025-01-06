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

    mth = zeros(1,Nfiles(l));       % Momentum thickness
    ruu = zeros(np(l),Nfiles(l));   % Reynolds stress: sqrt( rho u u)
    rvv = zeros(np(l),Nfiles(l));   % Reynolds stress: sqrt( rho v v)
    rww = zeros(np(l),Nfiles(l));   % Reynolds stress: sqrt( rho w w)
    ruv = zeros(np(l),Nfiles(l));   % Reynolds stress: sqrt(-rho u v)
    u_mean = zeros(np(l),Nfiles(l));% Mean streamwise velocity

    for q = 1:Nfiles(l)
        % Read data
        filename = strcat(mfc_dir(l),"/restart_data/lustre_",int2str(timesteps(l,q)),".dat");
        disp(filename);
        fileID = fopen(filename,'r');
        A = fread(fileID,'double');
        fclose(fileID);

        % Reassign density & velocity components
        rho     = zeros(  mp(l),np(l),pp(l));
        mom     = zeros(3,mp(l),np(l),pp(l));
        vel     = zeros(3,mp(l),np(l),pp(l));
        E       = zeros(  mp(l),np(l),pp(l));
        pres    = zeros(  mp(l),np(l),pp(l));
        alpha_b = zeros(  mp(l),np(l),pp(l));
        nbub    = zeros(  mp(l),np(l),pp(l));
        R       = zeros(nb,mp(l),np(l),pp(l));
        V       = zeros(nb,mp(l),np(l),pp(l));

        for k = 1:pp(l)
            for j = 1:np(l)
                for i = 1:mp(l)
                    % Continuity
                    rho(i,j,k)   = A(i+(j-1)*mp(l)+(k-1)*mp(l)*np(l));   % rho
                    % Momentum
                    mom(1,i,j,k) = A((momxb - 1)*r(l)+i+(j-1)*mp(l)+(k-1)*mp(l)*np(l)) * (u_b(l)/u_m(l));   % vel_x
                    mom(2,i,j,k) = A((momxb    )*r(l)+i+(j-1)*mp(l)+(k-1)*mp(l)*np(l)) * (u_b(l)/u_m(l));   % vel_y
                    if p > 0
                        mom(3,i,j,k) = A((momxe - 1)*r(l)+i+(j-1)*mp(l)+(k-1)*mp(l)*np(l)) * (u_b(l)/u_m(l));   % vel_z
                    end
                    vel(:,i,j,k) = mom(:,i,j,k)/rho(i,j,k);

                    if (bubbles == "T")
                        dyn_p = 0.5 * (mom(1,i,j,k)*vel(1,i,j,k) + mom(2,i,j,k)*vel(2,i,j,k) + mom(3,i,j,k)*vel(3,i,j,k));
                        % Energy
                        E(i,j,k) = A((E_idx - 1)*r(l)+i+(j-1)*mp(l)+(k-1)*mp(l)*np(l)) * (u_b(l)/u_m(l))^2;   % vel_z
    
                        % Void fraction
                        alpha(i,j,k) = A((advxe - 1)*r(l)+i+(j-1)*mp(l)+(k-1)*mp(l)*np(l));
                        % Number density
                        nbub(i,j,k) = A((n_idx - 1)*r(l)+i+(j-1)*mp(l)+(k-1)*mp(l)*np(l));
                        % nR, nV
                        for s = 1:nb
                            R(s,i,j,k) = A((bubxb - 1 + 2*(s - 1))*r(l)+i+(j-1)*mp(l)+(k-1)*mp(l)*np(l))/nbub(i,j,k);
                            V(s,i,j,k) = A((bubxb     + 2*(s - 1))*r(l)+i+(j-1)*mp(l)+(k-1)*mp(l)*np(l))/nbub(i,j,k);
                        end
                        % Void fraction from number density
                        R3 = 0;
                        for s = 1:nb
                            R3 = R3 + weight(s)*R(s,i,j,k)^3;
                        end
                        alpha(i,j,k) = (4*pi/3)*nbub(i,j,k)*R3;
                        % Pressure
                        pres(i,j,k) = (1/Gamma_w) * ((1/(1-alpha(i,j,k))) * (E(i,j,k) - dyn_p) - Pi_inf_w(l));
                    end
                end
            end
        end
        
        % Compute mean and fluctuating quantities
        rho_mean = mean(rho,[1 3]); % Mean density
        mom_mean = squeeze(mean(mom,[2 4])); % Mean momentum
        vel_mean = mom_mean./rho_mean; % Favre-averaged velocity

        % Compute velocity fluctuation
        rho_mean_rep = repmat(rho_mean,[mp(l) 1 pp(l)]);
        vel_mean_rep = permute(repmat(vel_mean,[1 1 mp(l) pp(l)]),[1 3 2 4]);
        rho_fluc = rho - rho_mean_rep;
        vel_fluc = vel - vel_mean_rep;

        % Compute vorticity thickness


        % Compute momentum thickness
        f = zeros(np(l),1); % integrand
        for j = 1:np(l)
            f(j) = rho_mean(j)*(1 - mom_mean(1,j)/rho_mean(j))*(1 + mom_mean(1,j)/rho_mean(j))/4;
        end
        mth(q) = dy(l)*trapz(f);

        % Normalized y-axis 
        y_norm = y(l,1:np(l))/mth(q);

        % Mean streamwise velocity
        u_mean(:,q) = vel_mean(1,:);

        if (Reynolds_stress == "T")
            % Compute Reynolds stress
            ruu(:,q) = squeeze(mean(rho_mean_rep.*squeeze(vel_fluc(1,:,:,:)).^2,[1 3]));
            rvv(:,q) = squeeze(mean(rho_mean_rep.*squeeze(vel_fluc(2,:,:,:)).^2,[1 3]));
            rww(:,q) = squeeze(mean(rho_mean_rep.*squeeze(vel_fluc(3,:,:,:)).^2,[1 3]));
            ruv(:,q) = squeeze(mean(rho_mean_rep.*squeeze(vel_fluc(1,:,:,:).*vel_fluc(2,:,:,:)),[1 3]));
        end

        if (TKE == "T")
 
            % Compute mean velocity gradient and velocity fluctuation gradient

            % Compute TKE terms - Transport, Production, Dissipation, Pressure-strain, mass flux coupling




        end

        if (bubbles == "T")
            if (ismember(timesteps(l,q),Nt_compare))
                flag = 1;
                for k = 1:pp(l)
                    for j = 1:np(l)
                        for i = 1:mp(l)
                            if ((abs(x(l,i) - bubbles_wrt_x(l)) < 0.5*dx(l)) && (abs(y(l,j) - bubbles_wrt_y(l)) < 0.5*dy(l)) && (abs(z(l,k) - bubbles_wrt_z(l)) < 0.5*dz(l)))
                                print_bubble_radius(strcat(p_bubble_radius_dir(l),int2str(timesteps(l,q)),".dat"), R(:,i,j,k)); disp('bubble_radius printed');

                                f_bubble_radius_ratio = figure("DefaultAxesFontSize", 18); plot_bubble_radius_ratio(f_bubble_radius_ratio, R(:,i,j,k));
                                saveas(f_bubble_radius_ratio, strcat(f_bubble_radius_ratio_dir(l),int2str(timesteps(l,q)),".png")); close(f_bubble_radius_ratio); disp('bubble_radius_ratio saved');

                                f_bubble_radius_pdf = figure("DefaultAxesFontSize", 18); plot_bubble_radius_pdf(f_bubble_radius_pdf, R(:,i,j,k));
                                saveas(f_bubble_radius_pdf, strcat(f_bubble_radius_pdf_dir(l),int2str(timesteps(l,q)),".png")); close(f_bubble_radius_pdf); disp('bubble_radius_pdf saved');
                                
                                flag = 0;
                                break;
                            end
                        end
                        if (flag == 0) 
                            break;
                        end
                    end
                    if (flag == 0) 
                        break;
                    end
                end
            end
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