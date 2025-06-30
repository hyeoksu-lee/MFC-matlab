function set_options(ncase)

    load variables/multiscale.mat;

    % Options
    read_raw_data       = "T";
    pres_stat           = "T";
    vorticity           = "F";
    Reynolds_stress     = "T";
    kolmogorov          = "F";
    tke_budget          = "F";
    energy_spectrum     = "F";
    avg_over_self_sim   = "F"; start_idx = 60;
    
    % Bubbles
    bubbles             = "F";
    bubbles_wrt_x       = 13404.64 * (x_b./x_m);
    bubbles_wrt_y       = -57.0410 * (x_b./x_m);
    bubbles_wrt_z       = 3365.420 * (x_b./x_m);

    save variables/options.mat;
end