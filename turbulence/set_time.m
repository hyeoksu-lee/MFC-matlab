function set_time(ncase)

    load variables/multiscale.mat;

    % Initialize variables
    Nt_beg = zeros(ncase,1);
    Nt_end = zeros(ncase,1);
    Nt_save = zeros(ncase,1);
    dt = zeros(ncase,1);

    % Time
    Nt_beg(:) = 0;
    Nt_end(:) = 1408*30;
    Nt_save(:) = 1408;
    dt(:) = 0.007102272727272727 .* (x_b./x_m) .* (u_m./u_b);

    % Save
    save variables/time.mat;
end
