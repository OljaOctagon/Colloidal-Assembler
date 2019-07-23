
module purge
module load gcc/7.2
module load boost/1.64.0
module load gsl/2.4 

#make a list of all files in directory
dir_arr=( `ls -d ./*/` )
counter=0

for i in "${dir_arr[@]}"
do
    echo $i >> batch_files.tmp
    counter=$(echo "$counter+1" | bc )
    if (( counter==16 ))
    then
        echo "new batch "
        counter=0
        sbatch submit_batch.slrm
        rm "batch_files.tmp"
    fi 

done
