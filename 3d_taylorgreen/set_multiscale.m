function set_multiscale(ncase)

    %% Air property
    rho_a = 1.2041;
    gamma_a = 1.4;
    pi_inf_a = 0;
    Gamma_a = 1/(gamma_a - 1);
    Pi_inf_a = gamma_a*pi_inf_a/(gamma_a - 1);

    %% Background flow property -- CAUTION
    gamma = gamma_a;
    pi_inf = pi_inf_a;
    Gamma = Gamma_a;
    Pi_inf = Pi_inf_a;

    %% Reference parameters
    P0 = 101325;
    rho0 = 1;
    C0 = sqrt(1.4 * P0);
    V0 = 0.1 * C0;
    L0 = 1;

    %% Reynolds number
    Re0 = 1600;

    save variables/multiscale.mat Gamma Pi_inf Re0 V0 P0 rho0 L0
end