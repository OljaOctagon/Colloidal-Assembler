import numpy as np 
import glob 
import pandas as pd 
import argparse 
from spanning_func import get_spanning
import gel_tools as gt 
import re
from scipy.optimize import curve_fit 
import pandas as pd 
import configparser
import matplotlib.pyplot as plt 


def plot_energy(energy_to_time, trend, level, fluct, is_converged, subdir_name):

    x = energy_to_time[:,0] 
    y = energy_to_time[:,1]

    fit_func = trend*x+ level 

    fig, ax = plt.subplots()

    plt.plot(x,y,lw=2, c='k')
    plt.plot(x,fit_func,c='r',lw=2, linestyle='dashed')
    plt.fill_between(x, fit_func+y, fit_func-y,facecolor='y', alpha=0.5)
    
    plt.savefig("energy_fit_{}.pdf".format(subdir_name))
    plt.close()

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


def get_energy_trend(energy_to_time):

    # get moving average 
    window_size = 6

    energy_ma = moving_average(energy_to_time[:,1],window_size)

    print("moving average")
    print(len(energy_ma))
    print(len(energy_to_time))

    # fit last n observations to linear function 
    def fit_func(x,k,d):
        return k*x + d

    last_nobs=20
    ydata=energy_ma[-20:]

    max_time = energy_to_time[-1,0]
    print("max time", max_time)
    freq=energy_to_time[1,0] - energy_to_time[0,0]
    print("freq", freq)


    xdata = np.linspace(max_time-freq*last_nobs,max_time,last_nobs)
    popt, pcov = curve_fit(fit_func,xdata,ydata)

    trend, level = popt

    fit_x = energy_to_time[window_size//2:-window_size//2,0]

    fluct = np.std(energy_to_time[-20:,1] - fit_func(energy_to_time[-20:,0],trend,level))

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
        pos = np.fromfile(pos_file)
        pos = np.reshape(pos, (-1,3))
        pos = pos[:,:2]

        box_file = "{}Box_{}.bin".format(dir_i,last_time)
        box = np.fromfile(box_file)

        yield ptype, phi, temperature, delta, last_time, pos, box


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
    

    gen = gen_dict[args.input]

    for ptype, phi, temperature, delta, last_time, pos, box in gen:
        new_results = {}
        # connections of last time 

        dir_name = "{}/{}_phi_{}_delta_{}_temp_{}".format(ptype,ptype,phi,delta,temperature)
        file_name = "{}/patch_network.dat".format(dir_name)

        subdir_name = "{}_phi_{}_delta_{}_temp_{}".format(ptype,phi,delta,temperature)

        connections = gt.read_bonds(file_name)[-1]   
        # calculate spanning 
        frac_largest, virtual_frac_largest = get_spanning(pos, box, connections)


        # Energy: trend and fluctuation estimate of last time points 
        energy_to_time = pd.read_csv("{}/Energy.dat".format(dir_name), delim_whitespace=True).values
        trend, level, fluct, is_converged = get_energy_trend(energy_to_time)

        plot_energy(energy_to_time, trend, level, fluct, is_converged,subdir_name)

        
        new_results['id'] = subdir_name 
        new_results['frac_largest'] = frac_largest 
        new_results['virtual_frac_largest'] = virtual_frac_largest
        new_results['energy_converged'] = is_converged 
        new_results['energy_fluctuation'] = fluct
        new_results['current_time'] = last_time 

        df = df.append(new_results, ignore_index=True)
        

    df.to_pickle("results_percolation.pickle")

