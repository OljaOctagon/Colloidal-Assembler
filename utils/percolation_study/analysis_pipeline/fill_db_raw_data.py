from psycopg2 import connect, extensions, sql
import psycopg2
import glob 
import configparser 
import re
import pickle
import numpy as np 

DB_NAME="test"
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

dirs = glob.glob("double*/")
print(dirs)

for dir_i in dirs: 

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

    mc_times = sorted(list(map(g,pos_files)))
    print("TIME",mc_times)

    for time_i in mc_times:


        pos_file = "{}positions_{}.bin".format(dir_i,time_i)
        with open(pos_file,'rb') as fh:
            pos = fh.read() 

        orient_file = "{}orientations_{}.bin".format(dir_i,time_i)
        with open(pos_file,'rb') as fh:
            orient = fh.read() 

        box_file = "{}Box_{}.bin".format(dir_i,time_i)
        with open(box_file,'rb') as fh:
            box = fh.read() 

        sql_statement="INSERT INTO data VALUES (%s,%s,%s,%s,%s,%s,%s,%s)"
        data_tuple = (ptype, delta, phi, temperature, time_i, pos, orient, box)
      
        #sql_statement="INSERT INTO data VALUES (%s,%s,%s,%s,%s)"
        #data_tuple = (ptype, delta, phi, temperature, time_i)

        cursor.execute(sql_statement, data_tuple)
        #cursor.execute('''UPDATE data SET pos = %s''', (pos + b'1',))
        #cursor.execute('''UPDATE data SET orient = %s''', (orient + b'1',))
        #cursor.execute('''UPDATE data SET box = %s''', (box + b'1',))


conn.close()

# TEST 
print("TEST TEST TEST TEST ")
print("test if binary data was read correctly")

expl_file="double_manta_asymm_1_phi_0.01_delta_0.2_temp_0.01/positions_20000000.bin" 

pos_from_buffer=np.fromfile(expl_file, dtype=float)
print(pos_from_buffer[0], len(pos_pos_from_buffer))


conn = connect(
dbname = DB_NAME,
user = "drcarina",
host = "localhost", 
password = "pwd"
)

sql="SELECT pos FROM data WHERE (phi,mctime) IN ((0.01, 20000000));"
cursor.execute(sql)
pos_from_db = pickle.loads(cursor.fetchall())
print(pos_from_db[0], len(pos_from_db))

