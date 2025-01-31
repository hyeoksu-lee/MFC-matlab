function print_vor_thickness(p_vor_thickness_dir, time, vth)

    fileID = fopen(p_vor_thickness_dir,'w');
    fprintf(fileID,'%12.8f %12.8f\r\n',[time; vth]);
    fclose(fileID);

end