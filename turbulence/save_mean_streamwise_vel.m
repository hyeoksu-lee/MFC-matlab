function save_mean_streamwise_vel(post_stat_dir, time, y_norm, u_mean)

    save(post_stat_dir+"/mean_streamwise_vel.mat","time","y_norm","u_mean");
    
end