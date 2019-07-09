for particle_type in double_manta double_mouse checkers
do
    for patch_type in symm_1 symm_2 Asymm_1 Asymm_2
    do
        cp calc_cluster_size.py $particle_type'/'$patch_type
        cd $particle_type/$patch_type
        rm psi_mean.json
        python calc_cluster_size.py -t $patch_type
        cp psi_mean.json ../../psi_mean_$particle_type"_"$patch_type.json
        cd ../..
    done
done
