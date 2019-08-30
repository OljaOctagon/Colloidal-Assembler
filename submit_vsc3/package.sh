
module purge
module load gcc/7.2 openmpi/3.0.1 python/2.7 boost/1.64.0
module load gsl/2.4 

#make a list of all files in directory
dir_arr=( `ls -d ./*/` )
counter=0
current_dir=$( echo $PWD )
echo "$PWD"
k=0
for i in "${dir_arr[@]}"
do
    echo $i >> $current_dir/batch_files_$k".tmp"
    counter=$(echo "$counter+1" | bc )
    if (( counter==20 ))
    then
        echo "new batch "
	counter=0
        sbatch submit_batch.slrm $current_dir/batch_files_$k".tmp" &
        wait %1
	k=$( echo "$k+1" | bc ) 
    fi 

done
