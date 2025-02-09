function print_min_pressure(p_min_pressure_dir, time, pres_min)

    fileID = fopen(p_min_pressure_dir,'w');
    fprintf(fileID,'%12.8f %12.8f\r\n',[time; pres_min]);
    fclose(fileID);

end