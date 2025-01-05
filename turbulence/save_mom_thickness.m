function save_mom_thickness(post_stat_dir, time, mth)

    save(post_stat_dir+"/momentum_thickness.mat","time","mth");
    
end