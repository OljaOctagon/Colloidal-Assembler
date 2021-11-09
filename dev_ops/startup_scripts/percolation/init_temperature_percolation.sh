mkdir init_runs

for d in ./double*/
   do
   name=${d#*/}
   name=${name%/}

   mv $name init_runs/
   echo "$d"
   ls init_runs/$name/RANDOM_STATE_R1_* > list.dat
   sed -i "s/.*_//" list.dat
   sed -i "s/\..*//" list.dat
   last_frame=$( awk 'BEGIN{max=0}{if ($1>max){max=$1}}END{print max}' < list.dat)

   for t in 0.01 0.04 0.05 0.07 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16
      do
      dirf=$name"_temp_"$t
      mkdir $dirf
      cp init_runs/$name"/positions_"$last_frame".bin" $dirf
      cp init_runs/$name"/orientations_"$last_frame".bin" $dirf
      cp init_runs/$name"/patch_energy_"$last_frame".bin" $dirf
      cp init_runs/$name"/Box_"$last_frame".bin" $dirf
      cp init_runs/$name"/RANDOM_STATE_R1_"$last_frame".bin" $dirf
      cp init_runs/$name"/RANDOM_STATE_R_"$last_frame".bin" $dirf
      cp init_runs/$name/McPoly $dirf
      cp init_runs/$name/para.ini $dirf
      sed -i "s/Temperature = 1/Temperature = $t/g" $dirf/para.ini
      done
done
~                   