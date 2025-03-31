function set_space(ncase)
    
    load variables/multiscale.mat;
    
    % Initialize space and time variables
    Lx_list = zeros(ncase,1); m_list = zeros(ncase,1);
    Ly_list = zeros(ncase,1); n_list = zeros(ncase,1);
    Lz_list = zeros(ncase,1); p_list = zeros(ncase,1);
     

    % x-direction
    Lx_list(:) = 160.0;
    m_list(:) = 511;
    stretch_x = "F";

    % y-direction
    Ly_list(:) = 160.0; 
    n_list(:) = 511; 
    stretch_y = "F"; a_y = 3; y_a = -0.3*Ly_list; y_b = 0.3*Ly_list; loops_y = 2;

    % z-direction
    Lz_list(:) = 80.0;
    p_list(:) = 255;
    stretch_z = "F";

    save variables/space.mat
end