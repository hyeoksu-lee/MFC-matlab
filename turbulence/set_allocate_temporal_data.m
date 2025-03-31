function set_allocate_temporal_data()

    load variables/global_parameters.mat
    load variables/options.mat

    %% Mean streamwise velocity
    u_mean = zeros(np,Nfiles); % Mean streamwise velocity

    %% Mixing layer thickness
    vth = zeros(1,Nfiles);       % Vorticity thickness
    mth = zeros(1,Nfiles);       % Momentum thickness

    %% y-normalized
    y_norm = zeros(np,Nfiles);

    %% Minimum pressure
    pres_min = zeros(1,Nfiles);

    %% Reynolds stress
    ruu = zeros(np,Nfiles);      % Reynolds stress: (rho u u)
    rvv = zeros(np,Nfiles);      % Reynolds stress: (rho v v)
    rww = zeros(np,Nfiles);      % Reynolds stress: (rho w w)
    ruv = zeros(np,Nfiles);      % Reynolds stress: (rho u v)

    %% Kolmogorov length scale
    eta_min = zeros(1,Nfiles);
    eta_max = zeros(1,Nfiles);

    %% Save variables
    save variables/temporal_data.mat
end