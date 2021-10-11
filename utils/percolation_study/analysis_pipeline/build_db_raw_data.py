#!/usr/bin/python3
# -*- coding: utf-8 -*-

# import the psycopg2 database adapter for PostgreSQL
from psycopg2 import connect, extensions, sql
import psycopg2

# declare a new PostgreSQL connection object
conn = connect(
dbname = "python_test",
user = "drcarina",
host = "localhost", 
password = "pwd"
)


# string for the new database name to be created
DB_NAME = "test"

# get the isolation leve for autocommit
autocommit = extensions.ISOLATION_LEVEL_AUTOCOMMIT
print ("ISOLATION_LEVEL_AUTOCOMMIT:", extensions.ISOLATION_LEVEL_AUTOCOMMIT)

"""
ISOLATION LEVELS for psycopg2
0 = READ UNCOMMITTED
1 = READ COMMITTED
2 = REPEATABLE READ
3 = SERIALIZABLE
4 = DEFAULT
"""

# set the isolation level for the connection's cursors
# will raise ActiveSqlTransaction exception otherwise
conn.set_isolation_level( autocommit )

# instantiate a cursor object from the connection
cursor = conn.cursor()

# use the execute() method to make a SQL request

cursor.execute("SELECT datname FROM pg_database;")
list_database = cursor.fetchall()

if (DB_NAME,) in list_database:
    print("'{}' Database already exists".format(DB_NAME))
else:
    print("'{}' Database does not exist, creating it...".format(DB_NAME))
    cursor.execute(sql.SQL("CREATE DATABASE {}").format(sql.Identifier( DB_NAME )))

# use the sql module instead to avoid SQL injection attacks

# close the cursor to avoid memory leaks
cursor.close()
# close the connection to avoid memory leaks
conn.close()


# reconnect and create table if it doesn't exist 

conn = connect(
dbname = DB_NAME,
user = "drcarina",
host = "localhost", 
password = "pwd"
)

conn.set_isolation_level( autocommit )

def table_exists(conn, TABLE_NAME):
    exists = False
    try:
        cursor = conn.cursor()
        cursor.execute("select exists(select relname from pg_class where relname='" + TABLE_NAME + "')")
        exists = cursor.fetchone()[0]
        cursor.close()
    except psycopg2.Error as e:
        print(e)
    return exists

TABLE_NAME="data"
exists=table_exists(conn, TABLE_NAME)
cursor = conn.cursor()

if exists:
    print("Data table already exists")

else:
    # create data table
    print("CREATE TABLE")
    #cursor.execute("CREATE TABLE data (ptype varchar, delta float, phi float, temperature float, mctime integer);")
    cursor.execute("CREATE TABLE data (ptype varchar, delta float, phi float, temperature float, mctime integer, pos bytea, orient bytea, box bytea, PRIMARY KEY (ptype, delta, phi, temperature, mctime));")
    #cursor.execute("CREATE TABLE data (ptype varchar, delta float, phi float, temperature float, mctime integer, pos bytea, orient bytea, box bytea);")


conn.close()
