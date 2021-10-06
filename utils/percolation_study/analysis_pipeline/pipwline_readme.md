
# in PSQL:
# Create use and database name 
psql
CREATE DATABASE python_test;
CREATE USER drcarina with encrypted password 'pwd_test';
GRANT ALL PRIVILEGES ON DATABASE python_test TO drcarina;
# Give user create priviledges
ALTER USER drcarina CREATEDB;

# python: create db for raw data  
python build_db_raw_data.py

