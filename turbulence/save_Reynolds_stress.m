function save_Reynolds_stress(post_stat_dir, time, y_norm, ruu, rvv, rww, ruv)

    save(post_stat_dir+"/Reynolds_stress.mat","time","y_norm","ruu","rvv","rww","ruv");
    
end