function set_user_inputs()

    %% Options
    fluid = "water";
    artificial_mach = true;
    Reynolds_stress = false;
    tke_budget = false;
    pres_stat = false;
    vorticity = false;
    energy_spectrum = false;
    get_liutex = true;
    growth_rate = false;

    %% Reynolds number
    Re = 320;   % Reynolds number
    Ca = 1.0;   % Cavitation number
    Ma_t = 0.1; % Target Mach number

    %% Grids
    Lx = 160.0; % Domain size in x-direction
    Ly = 160.0; % Domain size in y-direction
    Lz = 80.0;  % Domain size in z-direction

    m = 127;   % Number of grids in x-direction (0, ..., m)
    n = 127;   % Number of grids in y-direction (0, ..., n)
    p = 63;    % Number of grids in z-direction (0, ..., p)

    % Stretched grid in y-direction
    stretch_y = true; 
    a_y = 2; 
    y_a = -0.3*Ly; 
    y_b = 0.3*Ly; 
    loops_y = 2;

    %% Time
    tstep = 88; tskip = 10; dt = 0.056818181818181816;
    % Nt_beg = tstep*0;
    % Nt_end = tstep*40;
    % Nt_save = tstep*tskip;
    % timesteps = Nt_beg:Nt_save:Nt_end;
    timesteps = tstep*10;
    time = timesteps*(dt*tskip);
    Nfiles = length(timesteps);

    %% Configuration
    num_fluids = 1; % Number of fluids
    num_dims = 3;   % Number of dimensions

    %% Index
    contxb = 1;
    contxe = num_fluids;
    momxb  = contxe + 1;
    momxe  = contxe + num_dims;
    E_idx  = momxe + 1;
    advxb  = E_idx + 1;
    advxe  = E_idx + num_fluids;
    sys_size = advxe;

    %% Fluid properties
    % Air
    gamma_a  = 1.4;
    pi_inf_a = 0;
    % Water at 20degC
    rho_w = 1000;
    gamma_w  = 7.1;         % [1] specific heat ratio
    pi_inf_w = 306E+6;      % [N/m2] liquid stiffness
    pv       = 2.3388E+3;   % [N/m2] vapor pressure
    % Surrounding fluid
    if (fluid == "air")
        Gamma = 1/(gamma_a - 1);
        Pi_inf = gamma_a*pi_inf_a/(gamma_a - 1);
    elseif (fluid == "water")
        Gamma = 1/(gamma_w - 1);
        Pi_inf = gamma_w*pi_inf_w/(gamma_w - 1);
    end

    %% Flow configuration
    u_m  = 3.4343;                      % [m/s] upper stream velocity
    p0 = pv + Ca*(0.5*rho_w*u_m^2);     % [N/m2] Background pressure

    %% Artificial Mach number
    if (fluid == "water" && artificial_mach) 
        pi_fac = (rho_w * u_m^2 / (gamma_w * Ma_t^2) - p0) / pi_inf_w;
        Pi_inf = Pi_inf*pi_fac;
    end

    %% Non-dimensionalization
    Pi_inf = Pi_inf / (rho_w * u_m^2);
    pv = pv / (rho_w * u_m^2);

    %% Save user_inputs
    if ~exist("variables", 'dir')
        mkdir(strcat("variables"));
    end
    save variables/user_inputs.mat;

end
