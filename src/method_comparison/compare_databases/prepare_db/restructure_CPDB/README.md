# CellPhoneDB structure conversion

CellPhoneDB uses sql database, which is provided as .db, there are 6 tables, we extract these tables as csv. 

You can do it via sqlite tool, the below command extracts the gene_table and saves it as gene_table.csv

```sqlite3 -header -csv cellphone_pre.db "select * from gene_table;" > gene_table.csv```

Or you can use DB Browser for SQLite and, extract all the tables as csv through GUI. 


https://sqlitebrowser.org/


Detailed information can be found in the notebook. 