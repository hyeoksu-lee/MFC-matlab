function save_Reynolds_stress(post_stat_dir, time, vth, mth, ruu, rvv, rww, ruv)

    save(post_stat_dir+"/Reynolds_stress.mat","time","vth","mth","ruu","rvv","rww","ruv");
    
end