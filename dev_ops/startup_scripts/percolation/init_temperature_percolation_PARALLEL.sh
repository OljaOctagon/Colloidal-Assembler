startup_dir="dmo_as1_init_runs"
run_dir="parallel_runs_percolation" 
copy_dir="startup_settings"

n_sweeps_needed=2560000
n_freq=20000

flabel="double_mouse_asymm_1"

for d in $startup_dir/$flabel*/
   do
   echo "$d"
   ls $d/RANDOM_STATE_R1_* > list.dat
   sed -i "s/.*_//" list.dat
   sed -i "s/\..*//" list.dat
   last_frame=$( awk 'BEGIN{max=0}{if ($1>max){max=$1}}END{print max}' < list.dat)

   first_frame=$( "$last_frame-$n_sweeps_needed" | bc )
   current_frame=$first_frame

   for ptype in double_manta_asymm_1 double_mouse_asymm_1 double_mouse_symm_1 double_mouse_symm_2
   do 
      
      mkdir $run_dir/$ptype

      for t in 0.005 0.01 0.04 0.05 0.07 0.09 0.10 0.11 0.12 0.13 0.14 0.15 0.16
      
      do
      tail=value=${$d#*$flabel}
      dirf=$ptype$tail"_temp_"$t
      mkdir $dirf
      
      cp $d/"positions_"$current_frame".bin" $dirf
      cp $d"orientations_"$current_frame".bin" $dirf
      cp $d/"patch_energy_"$current_frame".bin" $dirf
      cp $d/"Box_"$current_frame".bin" $dirf
      cp $d/"RANDOM_STATE_R1_"$current_frame".bin" $dirf
      cp $d/"RANDOM_STATE_R_"$current_frame".bin" $dirf
      cp $d/"para.ini" $dirf
     
      # TODO Change more things here as ptypes have to be specified 
      cp $copy_dir/COMMITT.sh $dir_f
      cp $copy_dor/McPoly     $dir_f

      energy_level=-1

      sed -i "s/Temperature = 1/Temperature = $t/g" $dirf/para.ini
      sed -i "s/rhombus_type = $flabel/rhombus_type = $ptype/g" $dirf/para.ini
      sed -i "s/Energy_Level = 0/Energy_Level = $energy_level/g" $dirf/para.ini
   
      sed -i "s/Calculation_Frequency = 20000/Calculation_Frequency = 100000/g" $dirf/para.ini
      sed -i "s/Checkpoint_Frequency = 20000/Checkpoint_Frequency = 100000/g" $dirf/para.ini
      sed -i "s/Frame_Frequency = 2000000/Frame_Frequency = 20000000/g" $dirf/para.ini

      current_frame=$("$current_frame+$n_freq" | bc)

      done
   done
done
~                   