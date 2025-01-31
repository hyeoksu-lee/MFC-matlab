function print_kolmogorov_length(p_kolmogorov_length_dir, time, eta_max, eta_min)

    fileID = fopen(p_kolmogorov_length_dir,'w');
    fprintf(fileID,'%12.8f %12.8f %12.8f\r\n',[time; eta_max; eta_min]);
    fclose(fileID);

end