function set_time(ncase)

    load variables/multiscale.mat;

    % Initialize variables
    Nt_beg = zeros(ncase,1);
    Nt_end = zeros(ncase,1);
    Nt_save = zeros(ncase,1);
    Nfiles = zeros(ncase,1);

    dt = zeros(ncase,1);

    % Time
    Nt_beg(:) = 0;
    Nt_end(:) = 1112;
    Nt_save(:) = 8; 
    Nfiles(:)  = (Nt_end - Nt_beg)./Nt_save + 1;

    Nt_compare = [0:8:1112];

    dt(:) = 2.3620742079790262 * (x_b./x_m) .* (u_m./u_b);

    % 
    timesteps = zeros(ncase,max(Nfiles));
    time = zeros(ncase,max(Nfiles));
    for i = 1:ncase
        timesteps(i,1:Nfiles(i)) = Nt_beg(i):Nt_save(i):Nt_end(i);
        time(i,1:Nfiles(i)) = timesteps(i,1:Nfiles(i))*dt(i);
    end

    save variables/time.mat;
end