import numpy as np 
import glob 
import pandas as pd 
import argparse 
from spanning_func import get_spanning
import gel_tools as gt 

from scipy.optimize import curve_fit 
import pandas as pd 


def plot_energy(energy_to_time, trend, level, fluct, is_converged, dir_name):

	x = energy_to_time[:,0] 
	y = energy_to_time[:,1]

	fit_func = trend*x+ level 

	plt.plot(x,y,lw=2, c='k')
	plt.plot(x,fit_func,c='r',lw=2, linestyle='dashed')
	plt.fill_between(x, fit_func+y, fit_func-y,facecolor='y', alpha=0.5)
	
	plt.savefig("energy_fit_{}.pdf".format(dir_name))


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


def get_energy_trend(energy_to_time):

	# get moving average 
	window_size = 5
	energy_ma = moving_average(energy_to_time[:,1],window_size)

	# fit last n observations to linear function 
	def fit_func(k,d):
		return kx + d

	last_nobs=20
	ydata=energy_ma[:20]
	xdata = np.linspace(0,last_nobs,last_nobs)
	popt, pcov = curv_fit(fit_func,xdata, ydata)

	trend, level = popt
	fluct = np.std(energy_ma - fit_func(trend,level))

	T_trend = 0.01
	is_converged = False
	if trend < T_trend:
		is_converged = True 

	return trend, level, fluct, is_converged 


def generator_from_fsys(fsys_iterator):
	
	for dir_i in fsys_iterator:
		config = configparser.ConfigParser()
        config.read('{}para.ini'.format(dir_i))

        N = int(config['System']['Number_of_Particles'])
        phi = float(config['System']['Packing_Fraction'])
        temperature = float(config['System']['Temperature'])
        ptype = config['Rhombus']['rhombus_type']
        delta = config['Rhombus']['patch_delta']
        patch_size = config['Rhombus']['patch_size']

        pos_files =glob.glob('{}positions_*.bin'.format(dir_i))

        # get the last value from the string 
        g = lambda x: int(re.findall(r'\d+', x)[-1])

        mc_times = list(map(g,pos_files))

        last_time = np.max(mc_times)

        pos_file = "{}positions_{}.bin".format(dir_i,last_time)
        with open(pos_file,'rb') as fh:
            pos = fh.read() 


        orient_file = "{}orientations_{}.bin".format(dir_i,last_time)
        with open(pos_file,'rb') as fh:
            orient = fh.read() 

        box_file = "{}Box_{}.bin".format(dir_i,last_time)
        with open(box_file,'rb') as fh:
            box = fh.read() 

       
        # convert pos, orient box info to numpy arrays 

        yield ptype, phi, temperature, delta, last_time, pos, orient, box,


def generator_from_db():
	pass
    
	'''
    DB_NAME = "db_percol_raw_data"
    TABLE_NAME = "data"
    conn = connect(
    dbname = DB_NAME,
    user = "drcarina",
    host = "localhost", 
    password = "pwd"
    )

    # get the isolation leve for autocommit
    autocommit = extensions.ISOLATION_LEVEL_AUTOCOMMIT
    conn.set_isolation_level( autocommit )
    cursor = conn.cursor()

    #loop over different queries and yield data 

    sql="SELECT * FROM {} WHERE mctime IN {};".format(TABLE_NAME, last_time)
    cursor.execute(sql)
    row_data=cursor.fetchone()

    # WITH GROUP BY 


    yield ptype, phi, temperature, delta, last_time, pos, orient, box
 	''' 

def id_to_vars(id):
	pass

def vars_to_id(vars):
	pass


if __name__ == '__main__':
	
	# read data either through files system via glob or via db 
	parser = argparse.ArgumentParser()
	parser.add_argument('-input', type=str, choices=['fsys'])

	args = parser.parse_args()

	gen_fsys = generator_from_fsys(glob.glob("double*/double*/"))
	#gen_db = generator_from_db(some_iterator)

	gen_dict = {'fsys': gen_fsys}

 	columns = ['id', 'frac_largest', 'frac_largest_virtual', 'energy_converged', 'energy_fluctuation', 'current_time']
 	df = pd.DataFrame(columns=columns)

 
    for ptype, phi, temperature, delta, last_time, pos, orient, box, connections in gen(dict[args.input]):

    	new_results = {}
    	# connections of last time 

    	dir_name = "{}_phi_{}_delta_{}_temp_{}".format(ptype,phi,delta,temperature)
    	file_name = "{}/patch_network.dat".format(dir_name)

    	connections = gt.read_bonds(file_name)[-1]   
		# calculate spanning 
    	frac_largest, virtual_frac_largest = get_spanning(pos, box, connections)


    	# Energy: trend and fluctuation estimate of last time points 
    	energy_to_time = "{}/Energy.dat".format(dir_name)
    	trend, level, fluct, is_converged = get_energy_trend(energy_to_time)

    	plot_energy(energy_to_time, trend, level, fluct, is_converged)

    	
    	new_results['id'] = dir_name 
    	new_results['frac_largest'] = frac_largest 
    	new_results['virtual_frac_largest'] = virtual_frac_largest
    	new_results['energy_converged'] = is_converged 
    	new_results['energy_fluctuation'] = fluct
    	new_results['current_time'] = last_time 

    	df = df.append(new_results, ignore_index=True)
		

    df.to_pickle("results_percolation.pickle")

