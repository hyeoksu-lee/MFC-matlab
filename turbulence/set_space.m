function set_space(ncase)
    
    load variables/multiscale.mat;
    
    % Initialize space and time variables
    Lx = zeros(ncase,1);
    Ly = zeros(ncase,1);
    Lz = zeros(ncase,1);

    m = zeros(ncase,1);
    n = zeros(ncase,1);
    p = zeros(ncase,1);

    % x-direction
    Lx(:) = 59.0 * (x_b./x_m);
    m(:) = 399; 
    mp = m + 1; 
    m_max = max(m);
    mp_max = m_max + 1;
    dx = Lx./mp;

    % y-direction
    Ly = Lx * 3; 
    n(:) = 1199; 
    np = n + 1;
    n_max = max(n);
    np_max = n_max + 1;
    dy = Ly./np;

    % z-direction
    Lz = Lx;
    p = m;
    pp = p + 1;
    p_max = max(p);
    pp_max = p_max + 1;
    dz = Lz./pp;

    % grids
    x = zeros(ncase, mp_max); 
    y = zeros(ncase, np_max); 
    z = zeros(ncase, pp_max);
    r = mp.*np.*pp;

    for l = 1:ncase
        for i = 1:mp(l)
            x(l,i) = (i - 1)*dx(l) + 0.5*dx(l);
        end
        for j = 1:np(l)
            y(l,j) = (j - 1)*dy(l) + 0.5*dy(l) - 0.5*Ly(l);
        end
        for k = 1:pp(l)
            z(l,k) = (k - 1)*dz(l) + 0.5*dz(l);
        end
    end

    % Save variables
    save variables/space.mat
end