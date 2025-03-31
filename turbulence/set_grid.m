function set_grid()

    load variables/global_parameters.mat;

    %% Allocate memory
    x_cc = zeros(mp,1); x_cb = zeros(mpp,1);
    y_cc = zeros(np,1); y_cb = zeros(npp,1); 
    z_cc = zeros(pp,1); z_cb = zeros(ppp,1);

    %% x-dir
    for i = 1:mpp
        x_cb(i) = (i - 1)*Lx/mp;
    end
    if (stretch_x == "T")
        x_cb = x_cb / Lx;
        x_a = x_a / Lx;
        x_b = x_b / Lx;
        for j = 1:loops_x
            for i = 1:mpp
                x_cb(i) = x_cb(i)/a_x* ...
                                 (a_x + log(cosh(a_x*(x_cb(i) - x_a))) ...
                                      + log(cosh(a_x*(x_cb(i) - x_b))) ...
                                    - 2*log(cosh(a_x*(x_b - x_a)/2)));
            end
        end        
        x_cb = x_cb/(x_cb(mpp) - x_cb(1))*Lx;
    end
    dx = x_cb(2:mpp) - x_cb(1:mp);
    x_cc = (x_cb(2:mpp) + x_cb(1:mp))/2;

    %% y-dir
    for i = 1:npp
        y_cb(i) = (i - 1)*Ly/np - 0.5*Ly;
    end
    if (stretch_y == "T")
        y_cb = y_cb / Ly;
        y_a = y_a / Ly;
        y_b = y_b / Ly;
        for j = 1:loops_y
            for i = 1:npp
                y_cb(i) = y_cb(i)/a_y* ...
                                 (a_y + log(cosh(a_y*(y_cb(i) - y_a))) ...
                                      + log(cosh(a_y*(y_cb(i) - y_b))) ...
                                    - 2*log(cosh(a_y*(y_b - y_a)/2)));
            end
        end        
        y_cb = y_cb/(y_cb(npp) - y_cb(1))*Ly;
    end
    dy = y_cb(2:npp) - y_cb(1:np);
    y_cc = (y_cb(2:npp) + y_cb(1:np))/2;
    
    %% z-dir
    for i = 1:ppp
        z_cb(i) = (i - 1)*Lz/pp;
    end
    if (stretch_z == "T")
        z_cb = z_cb / Lz;
        z_a = z_a / Lz;
        z_b = z_b / Lz;
        for j = 1:loops_z
            for i = 1:ppp
                z_cb(i) = z_cb(i)/a_z* ...
                                 (a_z + log(cosh(a_z*(z_cb(i) - z_a))) ...
                                      + log(cosh(a_z*(z_cb(i) - z_b))) ...
                                    - 2*log(cosh(a_z*(z_b - z_a)/2)));
            end
        end        
        z_cb = z_cb/(z_cb(ppp) - z_cb(1))*Lz;
    end
    dz = z_cb(2:ppp) - z_cb(1:pp);
    z_cc = (z_cb(2:ppp) + z_cb(1:pp))/2;

    %% Save variables
    save variables/grid.mat x_cc y_cc z_cc x_cb y_cb z_cb dx dy dz
end