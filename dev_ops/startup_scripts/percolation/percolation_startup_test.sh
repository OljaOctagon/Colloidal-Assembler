sdir="/Users/ada/Documents/Code_Development_2020/rhombi/percolations_study/startup_all_particle_types_test"
newdir="/Users/ada/Documents/Code_Development_2020/rhombi/percolations_study/percolation_runs"

mkdir $newdir

for ptype in double_manta_asymm_1 double_mouse_symm_1 double_mouse_symm_2 double_mouse_asymm_1
do
    mkdir $newdir/$ptype
    for delta in  0.2 
    do
        ntime=1
        for temp in 0.001
        do
            for phi in  0.3  
	    do
	    ndir=$ptype"_temp_"$temp"_phi_"$phi"_delta_"$delta
            mkdir $newdir/$ptype/$ndir

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
	    
	    cd $newdir/$ptype/$ndir
	    bash COMMITT.sh & 
	    
    
    	    done
        done
    done
done


