function set_global_parameters(l)
    
    load variables/directory.mat;
    load variables/space.mat;
    load variables/time.mat;
    load variables/multiscale.mat;

    % Grid
    Lx = Lx_list(l); mp = m_list(l) + 1; mpp = mp + 1;
    Ly = Ly_list(l); np = n_list(l) + 1; npp = np + 1;
    Lz = Lz_list(l); pp = p_list(l) + 1; ppp = pp + 1;
    r = mp*np*pp;
    
    % Time steps
    Nfiles  = floor((Nt_end(l) - Nt_beg(l))/Nt_save(l)) + 1;
    timesteps = Nt_beg(l):Nt_save(l):Nt_end(l);
    time = timesteps*dt(l);

    % Post stat saving directory
    post_stat_dir = mfc_dir(l) + "/post_stat";

    % Reynolds number
    Re = Re0;

    save variables/global_parameters.mat
end
