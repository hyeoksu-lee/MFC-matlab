function set_index()

    num_fluids = 1;
    num_dims = 3;
    bubbles = "F";
    adv_n = "F";

    % Index
    contxb = 1;
    contxe = num_fluids;
    momxb  = contxe + 1;
    momxe  = contxe + num_dims;
    E_idx  = momxe + 1;
    advxb  = E_idx + 1;
    advxe  = E_idx + num_fluids;
    sys_size = advxe;

    save variables/index.mat;
end