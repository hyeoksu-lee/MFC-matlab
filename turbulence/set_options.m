function set_options(ncase)

    load variables/multiscale.mat;

    bubbles_wrt_x = zeros(ncase,1);
    bubbles_wrt_y = zeros(ncase,1);
    bubbles_wrt_z = zeros(ncase,1);

    % Options
    mom_thickness       = "F";
    mean_streamwise_vel = "F";
    Reynolds_stress     = "F";
    min_pressure        = "T";
    vortex              = "T";

    % Bubbles
    bubbles             = "F";
    bubbles_wrt_x       = 13404.64 * (x_b./x_m);
    bubbles_wrt_y       = -57.0410 * (x_b./x_m);
    bubbles_wrt_z       = 3365.420 * (x_b./x_m);

    save variables/options.mat;
end