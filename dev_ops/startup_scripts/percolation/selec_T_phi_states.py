import os 
import argparse
import pandas as pd 

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_agrument(
		"-input", 
		description="Pickle file that contains the largest virtual cluster for all systems. Naming scheme: results_percolation_<run_id>.pickle")

    args = parser.parse_args()
    df = pd.read_pickle(args.input)

    df = df.sort_values(['ptype', 'delta', 'phi','temperature'])


    def filter_function(T1,phi1,cluster_size1,cluster_size2):
    	low_density_cutoff = 0.227
    	is_percol_candidate = False 

    	if (phi1<=low_density_cutoff) and (T1<0.11):
    		is_percol_candidate = True

    	if (phi>low_density_cutoff):
    		if (cluster_size1 > 0.1) and (cluster_size1<0.5):
    			is_percol_candidate=True 

    		else:
    			if (cluster_size2 > 0.1) and (custer_size2 < 0.5):
    				is_percol_candidate=True 
    	return is_percol_candidate

   df['is_percol_candidate'] = df.apply(
   	filter_function(df.T,df.phi,df.frac_largest_virtual.shift(1)))

   run_id = df.id.values
   is_percol_candidate = df.is_percol_candidate.values 

   for r_i, percol_i in zip(run_id,is_percol_candidate):
       if percol_i == False:
           print("Going to remove {}".format(r_i))
           #os.rmdir("./{}".format(r_i))

  