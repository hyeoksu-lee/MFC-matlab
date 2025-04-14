function set_time(ncase)

    load variables/multiscale.mat;

    % Initialize variables
    Nt_beg = zeros(ncase,1);
    Nt_end = zeros(ncase,1);
    Nt_save = zeros(ncase,1);
    dt = zeros(ncase,1);

    % Time
    Nt_beg(:) = 0;
    Nt_end(:) = 1387*24;
    Nt_save(:) = 1387;
    dt(:) = 0.14911816459390956 .* (x_b./x_m) .* (u_m./u_b);

    % Save
    save variables/time.mat;
end
