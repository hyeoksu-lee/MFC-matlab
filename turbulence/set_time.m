function set_time(ncase)

    load variables/multiscale.mat;

    % Initialize variables
    Nt_beg = zeros(ncase,1);
    Nt_end = zeros(ncase,1);
    Nt_save = zeros(ncase,1);
    dt = zeros(ncase,1);

    % Time
    Nt_beg(:) = 496*25;
    Nt_end(:) = 496*40;
    Nt_save(:) = 496;
    dt(:) = 0.4175308608629467 .* (x_b./x_m) .* (u_m./u_b);

    % Save
    save variables/time.mat;
end
