function print_bubble_radius(p_bubble_radius_dir, radius)

    load variables/poly_bubbles.mat;

    fileID = fopen(p_bubble_radius_dir,'w');
    A = [R0 radius weight];
    fprintf(fileID,'%12.8e %12.8e %12.8e\r\n',A');
    fclose(fileID);

end