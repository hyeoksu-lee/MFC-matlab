function set_multiscale(ncase)

    x_m = zeros(ncase,1); u_m = zeros(ncase,1);
    x_b = zeros(ncase,1); u_b = zeros(ncase,1);
    pres0 = zeros(ncase,1);

    % Fluid property
    rho_w = 1000;               % [kg/m3] density of water
    gamma_w = 7.1;              % [1] specific heat ratio
    pi_inf_w = 306E+06 * 0.00208675338389119;         % [N/m2] liquid stiffness
    Gamma_w = 1/(gamma_w - 1);
    Pi_inf_w = gamma_w*pi_inf_w/(gamma_w - 1);

    rho_a = 1.2041;
    gamma_a = 1.4;
    pi_inf_a = 0;
    Gamma_a = 1/(gamma_a - 1);
    Pi_inf_a = gamma_a*pi_inf_a/(gamma_a - 1);

    % Mixing layer scale
    x_m(:) = 0.002475;          % [m] initial vorticity thickness
    u_m(:) = 6.8686;            % [m/s] upper-lower stream velocity difference

    pres0(:) = 25927.63298;      % [N/m2] Initial pressure corresponding to Ca = 1
    Pi_inf_w = Pi_inf_w ./ pres0;

    % Sub-grid bubble scale
    x_b(:) = 10E-6;             % [m] Initial bubble equilibrium radius
    u_b(:) = sqrt(pres0./rho_w); % [m/s] Reference velocity scale of bubble

    % For single phase flow
    x_m(:) = 1;
    u_m(:) = 1;
    x_b(:) = 1;
    u_b(:) = 1;
    Re0 = 160;
    
    save variables/multiscale.mat;
end