
data_location[0]=/global/lv70753/karner4/patchy_rhombi/FINAL_RUNS_17_OKT
data_location[1]=/global/lv70753/karner4/patchy_rhombi/FINAL_RUNS_17

keytype_5kbt[0]=batterie_1_5.2
keytype_5kbt[1]=batterie_2_5.2

keytype_kbt[0]=lege_batterie1
keytype_kbt[1]=lege_batterie2

############################

psi_data=psi_data
mkdir $psi_data

for particle_type in double_manta double_mouse checkers
do
    mkdir $psi_data/$particle_type
    mkdir $psi_data/$particle_type/symm_1
    mkdir $psi_data/$particle_type/symm_2
    mkdir $psi_data/$particle_type/Asymm_1
    mkdir $psi_data/$particle_type/Asymm_2
    for j in 0 1
    do
        k=$( echo "$j+1" | bc)
        cp ${data_location[0]}/$particle_type/${keytype_5kbt[$j]}/analysis/psi_op_mu_0.25*2symm* $psi_data/$particle_type/symm_$k
        cp ${data_location[0]}/$particle_type/${keytype_5kbt[$j]}/analysis/psi_op_mu_0.25*2Asymm* $psi_data/$particle_type/Asymm_$k

        cp ${data_location[1]}/$particle_type/${keytype_kbt[$j]}/analysis/psi_op_mu_0.25*2symm* $psi_data/$particle_type/symm_$k
        cp ${data_location[1]}/$particle_type/${keytype_kbt[$j]}/analysis/psi_op_mu_0.25*2Asymm* $psi_data/$particle_type/Asymm_$k

    done
done
