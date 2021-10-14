from psycopg2 import connect, extensions, sql
import psycopg2
import numpy as np 
import pickle 
import unittest

class Test_DB_values(unittest.TestCase):

    def test_reading_bytea(self):

        expl_file="double_manta_asymm_1_phi_0.01_delta_0.2_temp_0.01/positions_20000000.bin"
        pos_from_buffer=np.fromfile(expl_file, dtype=float)
       
        ###############################

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

        sql="SELECT pos FROM data WHERE (phi,mctime) IN ((0.01, 20000000));"
        cursor.execute(sql)
        pos_from_db = pickle.loads(cursor.fetchall())
       
        self.assertEqual(pos_from_db[0],pos_from_buffer[0])
        self.assertEqual(len(pos_from_db),len(pos_from_buffer))


if __name__ == '__main__':
    unittest.main()
