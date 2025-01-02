function print_mom_thickness(p_mom_thickness_dir, time, mth)

    fileID = fopen(p_mom_thickness_dir,'w');
    fprintf(fileID,'%12.8f %12.8f\r\n',[time; mth]);
    fclose(fileID);

end