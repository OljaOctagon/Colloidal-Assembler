from psycopg2 import connect, extensions, sql
import psycopg2
import glob 
import configparser 
import re

if __name__ == '__main__':
    
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

    dirs = glob.glob("double*/double*")
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
          
            cursor.execute(sql_statement, data_tuple)
       

    conn.close()

    