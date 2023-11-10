# load hgmd_phenbase data dump
# mysql create database
mysql> CREATE DATABASE hgmd_phenbase;
mysql> exit;

# load data dump
mysql -u USER -p hgmd_phenbase < hgmd_phenbase-2022.2.dump.sql

# query acc_num to HPO ID mapping (hgmd_phen.sql in GitHub repo under the same path as this one)
mysql -u USER -p hgmd_phenbase < hgmd_phen.sql

# set the output path in hgmd_phen.sql
# query result is saved at "/var/lib/mysql-files/HGMD_phen.tsv"