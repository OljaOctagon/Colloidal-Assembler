import os
import re
import numpy as np
import pandas as pd
from configparser import SafeConfigParser
parser = SafeConfigParser()
directory_list = [
'mu_0.25Energy_-5.2Asymm_patchpos_0.2_1',
'mu_0.25Energy_-5.2Asymm_patchpos_0.3_1'
'mu_0.25Energy_-5.2Asymm_patchpos_0.4_1n']

#directory_list = os.listdir('.')

for dir in directory_list: 

    files = os.listdir(dir)
    checkpoint_list = [ re.findall('\d+', file)[0] for file in files if file.startswith("pos") and file.endswith("bin") ]
    max_cp = max(list(map(int,checkpoint_list)))
    #os.system('rm -r '+dir+'/orderparameter')
    os.system('mkdir '+dir+'/orderparameter')

    os.system('cp '+dir+'/positions_'+str(max_cp)+'.bin '+dir+'/orderparameter')
    os.system('cp '+dir+'/orientations_'+str(max_cp)+'.bin '+dir+'/orderparameter')
    os.system('cp '+dir+'/Box_'+str(max_cp)+'.bin '+dir+'/orderparameter')
    os.system('cp '+dir+'/RANDOM_STATE_R1_'+str(max_cp)+'.bin '+dir+'/orderparameter')
    os.system('cp '+dir+'/RANDOM_STATE_R_'+str(max_cp)+'.bin '+dir+'/orderparameter')

    os.system('cp McPoly '+dir+'/orderparameter')
    os.system('cp para.ini '+dir+'/orderparameter')
    
    mu, energy, patch_delta, run = re.findall(r"[-+]?\d*\.\d+|\d+", dir)
    
    parser.read(dir+'/orderparameter/para.ini')
    parser.set('Rhombus', 'patch_delta', patch_delta)
    parser.set('Rhombus', 'Energy_Level', energy)
    parser.set('System', 'mu', mu)

    os.system('cd '+dir+'/orderparameter; ./McPoly -f '
        +str(max_cp)+' '+str(max_cp+1000))

    print(energy)
    arr = pd.read_csv(dir+'/orderparameter/psi_op.dat', 
        delim_whitespace=True).values[:,1]
    mean = np.mean(arr)
    std = np.std(arr)

    with open("double_manta_asymm_2_psi.dat", 'a+') as f:
        f.write(str(patch_delta)+' '+str(mean)+' '+str(std))


