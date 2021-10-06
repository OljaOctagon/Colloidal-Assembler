from psycopg2 import connect, extensions, sql
import psycopg2
import glob 


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

files_i = "double*/*.bin"


# order files in right order 
# split and find out identifier 




files = glob.glob(ndir)

        mean_psi = []
        std_psi = []

        mean_N = []
        std_N = []

        for file_i in files:

            phi_i = file_i.split("_")[2]
            temp_i = file_i.split("_")[6]  







sql_statement="INSERT INTO data VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)"
data_tuple = 


cursor.execute(sql_statement, data_tuple)

cursor.execute("ALTER TABLE data ADD ROW ")

cursor.execute("CREATE TABLE data (ptype varchar, delta float, phi float, temperature float,  prun integer, mctime integer, pos bytea, orient bytea, box bytea, PRIMARY KEY (ptype, delta, phi, temperature, prun, mctime));")

