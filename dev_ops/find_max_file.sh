for f in $(pwd)/*  
do 
 if [ -d $f ] 
   then
	ls $f/RANDOM_STATE_R1_* > list.dat
	sed -i "s/.*_//" list.dat	
 	sed -i "s/\..*//" list.dat
	b=$( awk 'BEGIN{max=0}{if ($1>max){max=$1}}END{print max}' < list.dat)
	a=$( echo "$b+3000000" | bc ) 
	cp COMMITT.sh $f
	cp McPoly $f
	sed -i "s/startpoint/$b/"  $f/COMMITT.sh
	sed -i "s/endpoint/$a/" $f/COMMITT.sh 

 fi 
done
