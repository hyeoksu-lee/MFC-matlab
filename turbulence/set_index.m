function set_index()

    load variables/poly_bubbles.mat;

    num_fluids = 1;
    num_dims = 3;
    bubbles = "T";
    adv_n = "T";

    % Index
    contxb = 1;
    contxe = num_fluids;
    momxb  = contxe + 1;
    momxe  = contxe + num_dims;
    E_idx  = momxe + 1;
    advxb  = E_idx + 1;
    advxe  = E_idx + num_fluids;
    sys_size = advxe;
    if (bubbles == "T")
        bubxb  = sys_size + 1;
        bubxe  = sys_size + 2*nb;
        sys_size = bubxe;
        if (adv_n == "T")
            n_idx  = sys_size + 1;
            sys_size = n_idx;
        end
    end

    save variables/index.mat;
end