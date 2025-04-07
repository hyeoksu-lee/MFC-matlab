function set_time(ncase)

    load variables/multiscale.mat;

    % Initialize variables
    Nt_beg = zeros(ncase,1);
    Nt_end = zeros(ncase,1);
    Nt_save = zeros(ncase,1);
    dt = zeros(ncase,1);

    % Time
    Nt_beg(:) = 667*0;
    Nt_end(:) = 667*80;
    Nt_save(:) = 667;
    dt(:) = 0.15509874887191297 .* (x_b./x_m) .* (u_m./u_b);

    % Save
    save variables/time.mat;
end
