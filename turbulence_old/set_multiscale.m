function set_multiscale(ncase)

    x_m = zeros(ncase,1); u_m = zeros(ncase,1);
    x_b = zeros(ncase,1); u_b = zeros(ncase,1);
    pres0 = zeros(ncase,1);

    %% Water property
    rho_w = 1000;               % [kg/m3] density of water
    gamma_w = 7.1;              % [1] specific heat ratio
    pi_inf_w = 306E+06;         % [N/m2] liquid stiffness
    Gamma_w = 1/(gamma_w - 1);
    Pi_inf_w = gamma_w*pi_inf_w/(gamma_w - 1);
    pv = 2338.8;

    %% Air property
    rho_a = 1.2041;
    gamma_a = 1.4;
    pi_inf_a = 0;
    Gamma_a = 1/(gamma_a - 1);
    Pi_inf_a = gamma_a*pi_inf_a/(gamma_a - 1);

    %% Nondimensional parameters
    Re0 = 50;
    Ca0 = 1;
    Ma_t = 0.1; % Target Mach number
    
    %% Mixing layer scale
    x_m = 0.002475/2;        % [m] half of initial vorticity thickness
    u_m = 3.4343;          % [m/s] velocity difference
    pres0 = pv + Ca0*(0.5*rho_w*u_m^2); % [N/m2] Initial pressure

    %% Sub-grid bubble scale
    x_b = 50E-6;             % [m] Initial bubble equilibrium radius
    u_b = sqrt(pres0/rho_w); % [m/s] Reference velocity scale of bubble

    %% Artificial Mach number
    pi_fac = (rho_w * u_m^2 / (gamma_w * Ma_t^2) - pres0) / pi_inf_w;

    %% Background flow property -- CAUTION
    Gamma = Gamma_w;
    Pi_inf = Pi_inf_w*pi_fac;
    Pi_inf = Pi_inf / (rho_w*u_m^2);
    pv = pv / (rho_w*u_m^2);
    
    %% For single phase flow
    % x_m = 1;
    % u_m = 1;
    % x_b = 1;
    % u_b = 1;

    save variables/multiscale.mat Gamma Pi_inf pv Re0 x_m u_m x_b u_b;
end