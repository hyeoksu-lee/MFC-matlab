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
% Set options
set_options(ncase); load variables/options.mat;
% Set index
set_index(); load variables/index.mat;
toc

%% POST-PROCESS
if (read_raw_data == "T")
    for l = 1:ncase
        disp("setup each case")
        tic
        set_global_parameters(l); load variables/global_parameters.mat;

        set_grid(); load variables/grid.mat;

        Ek = zeros(Nfiles,1);
        diss1 = zeros(Nfiles,1);
        diss2 = zeros(Nfiles,1);

        toc

        % Parallized loop over timesteps
        for q = 1:Nfiles
            
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

            % Total kinetic energy
            Ek(q) = mean(0.5*qp(1,:,:,:).*(qp(2,:,:,:).^2 + qp(3,:,:,:).^2 + qp(4,:,:,:).^2),[2 3 4]);

            % Compute derivatives
            disp("vel_derivatives")
            tic
            dvel_ds = f_compute_vel_derivatives(qp(momxb:momxe,:,:,:));
            toc

            % Compute vorticity
            disp("vorticity")
            tic
            omega = f_compute_vorticity(dvel_ds);
            enstropy = qp(1,:,:,:).*( omega(1,:,:,:).^2 + omega(2,:,:,:).^2 + omega(3,:,:,:).^2)/2;
            diss2(q) = (2/Re)*mean(enstropy,[2 3 4]);
            toc

        end

        diss1(1) = -(Ek(2) - Ek(1))/(time(2) - time(1));
        diss1(end) = -(Ek(end) - Ek(end-1))/(time(end) - time(end-1));
        diss1(2:end-1) = -(Ek(3:end) - Ek(1:end-2))./(time(3:end) - time(1:end-2))';

        % SAVE DATA AS .MAT FILES
        disp("save data")
        tic
        save(post_stat_dir+"/dissipation_rate.mat","time","diss1","diss2");
        toc
    end
end

if (plot_compare == "T")
    f1 = figure("DefaultAxesFontSize",18);
    set(f1,"Position",[100 100 1200 600]);
    A = readmatrix("data/conny.dat");
    time_ref = A(:,1); diss_ref = A(:,2);
    subplot(1,2,1);
    plot(time_ref,diss_ref,'ko','LineWidth',2); hold on; grid on;
    axis([0 20 0 0.018]);
    xlabel('$t$','interpreter','latex');
    ylabel('$\epsilon$','interpreter','latex');

    subplot(1,2,2);
    plot(time_ref,diss_ref,'ko','LineWidth',2); hold on; grid on;
    axis([0 20 0 0.018]);
    xlabel('$t$','interpreter','latex');
    ylabel('$\epsilon$','interpreter','latex');

    for l = 1:ncase

        set_global_parameters(l); load variables/global_parameters.mat;
        load(post_stat_dir+"/dissipation_rate.mat");

        subplot(1,2,1);
        plot(time,diss1,'LineWidth',2); hold on;

        subplot(1,2,2);
        plot(time,diss2,'LineWidth',2); hold on;

    end
    legend("Reference","WENO3M, N=127", "WENO5M, N=127", "WCNS6LD, N=127");
    saveas(f1,output_dir);
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
    qc(momxb:momxe,:,:,:) = qc(momxb:momxe,:,:,:) / (rho0*V0);  % Momentum
    qc(E_idx,:,:,:) = qc(E_idx,:,:,:) / (P0);    % Energy

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
    qp_mean = zeros(sys_size,1);

    % mean rho
    qp_mean(1) = squeeze(mean(qp(1,:,:,:),[2 3 4]));

    % favre-averaged velocity
    for i = momxb:momxe
        qp_mean(i) = squeeze(mean(qp(1,:,:,:).*qp(i,:,:,:),[2 3 4])./mean(qp(1,:,:,:),[2 3 4]));
    end

    % mean pressure
    qp_mean(E_idx) = squeeze(mean(qp(E_idx,:,:,:),[2 3 4]));

    % fluctuation
    qp_fluc = qp - repmat(qp_mean, [1 mp np pp]);
end

% f_compute_vel_derivatives
function [dvel_ds] = f_compute_vel_derivatives(vel)

    load variables/global_parameters.mat;
    load variables/grid.mat;
    
    dvel_ds = zeros(3,3,mp,np,pp);

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
