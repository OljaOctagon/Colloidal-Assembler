from psycopg2 import connect, extensions, sql
import psycopg2
import glob 
import configparser 
import re

DB_NAME="raw_data_rhombi_percolation"
TABLE_NAME = "data"

conn = connect(
dbname = DB_NAME,
user = "drcarina",
host = "localhost", 
password = "pwd_test"
)

conn.set_isolation_level( autocommit )
cursor = conn.cursor()

dirs = glob.glob("double*")
print(dirs)

for dir_i in dirs: 

    config = configparser.ConfigParser()
    config.read('{}/para.ini'.format(dir_i))

    N = int(config['System']['Number_of_Particles']
    phi = float(config['System']['Packing_Fraction'])
    temperature = float(config['System']['Temperature'])
    ptype = config['Rhombus']['rhombus_type']
    delta = config['Rhombus']['patch_delta']
    patch_size = config['Rhombus']['patch_size']

    pos_file =glob.glob('{}/positions_*.bin')
    orient_file = glob.glob('{}/orientations_*.bin')
    box_file = glob.glob("{}/Box_*.bin")

    # get the last value from the string 
    g = lambda x: int(re.findall(r'\d+', x)[-1])

    mc_times = sorted(list(map(g,pos_file)))

    for time_i in mc_times:
        with open(pos_file,'rb') as fh:
            pos_bin = fh.read() 

        orient_file = "{}/orientations_{}.bin".format(dir_i,time_i)
        with open(pos_file,'rb') as fh:
            orient_bin = fh.read() 

        box_file = "{}/Box_{}.bin".format(dir_i,time_i)
        with open(box_file,'rb') as fh:
            box_bin = fh.read() 

       
        sql_statement="INSERT INTO data VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)"
        data_tuple = 

        cursor.execute(sql_statement, data_tuple)
        c.execute("INSERT INTO players(player_name) VALUES(%(name)s)", {"name":name})
        