function set_space(ncase)
    
    load variables/multiscale.mat;
    
    % Initialize space and time variables
    Lx_list = zeros(ncase,1); m_list = zeros(ncase,1);
    Ly_list = zeros(ncase,1); n_list = zeros(ncase,1);
    Lz_list = zeros(ncase,1); p_list = zeros(ncase,1);
     

    % x-direction
    Lx_list(:) = 2*pi; %160.0;
    m_list(:) = 127;
    stretch_x = "F";

    % y-direction
    Ly_list = Lx_list; %160.0;
    n_list = m_list; 
    stretch_y = "F";

    % z-direction
    Lz_list = Lx_list; %80.0;
    p_list = m_list;
    stretch_z = "F";

    save variables/space.mat
end