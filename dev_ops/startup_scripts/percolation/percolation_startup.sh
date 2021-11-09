sdir="/Users/ada/Documents/Code_Development_2020/rhombi/percolations_study/startup_all_particle_types"
newdir="/Users/ada/Documents/Code_Development_2020/rhombi/percolations_study/percolation_runs"

#sdir="/home/lv71286/karner_ipc/rectangles/rectangles_mono_startup"
#newdir="/home/lv71286/karner_ipc/rectangles/rectangles_mono"
mkdir $newdir

for ptype in double_manta_asymm_1 double_mouse_symm_1 double_mouse_symm_2 double_mouse_asymm_1
do
    mkdir $newdir/$ptype
    for delta in 0.1 0.2 0.3 0.4 0.5
    do
        ntime=1
        for temp in 0.001 0.01 0.02 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2
        do
            for phi in 0.005 0.1 0.2 0.3 0.4 0.5 
	    do
	    ndir=$ptype"_temp_"$temp"_phi_"$phi"_delta_"$delta
            mkdir $newdir/$ptype/$ndir

            #cp $sdir/$ptype/positions_$ntime".bin"        $newdir/$ptype/$ndir
            #cp $sdir/$ptype/orientations_$ntime".bin"     $newdir/$ptype/$ndir
            #cp $sdir/$ptype/Box_$ntime".bin"              $newdir/$ptype/$ndir
            #cp $sdir/$ptype/RANDOM_STATE_R_$ntime".bin"   $newdir/$ptype/$ndir
            #cp $sdir/$ptype/RANDOM_STATE_R1_$ntime".bin"  $newdir/$ptype/$ndir
            #cp $sdir/$ptype/patch_energy_$ntime".bin"     $newdir/$ptype/$ndir

            cp $sdir/COMMITT.sh $newdir/$ptype/$ndir
            cp $sdir/para.ini   $newdir/$ptype/$ndir
            cp $sdir/McPoly     $newdir/$ptype/$ndir

            gsed -i "s/ptype/$ptype/g"   $newdir/$ptype/$ndir/para.ini
            gsed -i "s/current_tmp/$temp/g"   $newdir/$ptype/$ndir/para.ini
            gsed -i "s/current_phi/$phi/g"   $newdir/$ptype/$ndir/para.ini
            gsed -i "s/p_value/$delta/g"   $newdir/$ptype/$ndir/para.ini

            a=$ntime
            b=$(echo "$a+10000000" | bc)

            gsed -i "s/startpoint/$a/"  $newdir/$ptype/$ndir/COMMITT.sh
            gsed -i "s/endpoint/$b/"    $newdir/$ptype/$ndir/COMMITT.sh
	    done
        done
    done
done


