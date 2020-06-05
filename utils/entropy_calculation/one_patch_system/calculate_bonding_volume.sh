for radius in 0.05 0.06 0.07 0.08 0.09 0.1 0.11 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.2 
do 
  for delta in 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8  
  do
	  fdir="bonding_volume_delta_"$delta"_radius_"$radius 
	  mkdir $fdir
	  cp McPoly $fdir
	  cp para.ini $fdir
	  cp make_parallel_rhombi.py $fdir 
          cp bonding_volume.py $fdir

	  gsed -i "s/pdelta/$delta/g" $fdir/para.ini 
	  gsed -i "s/pradius/$radius/g" $fdir/para.ini
	  
	  cd $fdir 
	  rand1=$RANDOM 
	  rand2=$RANDOM 
	  ./McPoly -n 0 1 $rand1 $rand2
	  rm patch_state.dat 
	  python make_parallell_rhombi.py -delta $delta
	  ./McPoly -f 1 100 $rand1 $rand2 
	  python bonding_volume.py -delta $delta -radius $radius
	  
	  cd ..  
	done 
done
