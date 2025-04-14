function set_time(ncase)

    load variables/multiscale.mat;

    % Initialize variables
    Nt_beg = zeros(ncase,1);
    Nt_end = zeros(ncase,1);
    Nt_save = zeros(ncase,1);
    dt = zeros(ncase,1);

    % Time
    % Nt_beg(:) = 0;
    % Nt_end(:) = 4000;%[4000 4000 4000];
    % Nt_save(:) = 40;
    % dt(:) = 0.00013033095413020053 * (V0/L0);

    Nt_beg(:) = 0;
    Nt_end(:) = [8100 8100 7209];
    Nt_save(:) = 81;
    dt(:) = 6.516547706510026e-05 * (V0/L0);

    % Save
    save variables/time.mat;
end
