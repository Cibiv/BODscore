#!/bin/bash
# script to convert `SQLite3` database produced by `vshape` to a `csv` file

dbfile="test.db"
table="testsample__coverage"
outputto="testtable.tab"
sep="\t"
fields="contig, pos, score"

echo "FOLLOWING TABLES ARE AVAILABLE:"

sqlite3 $dbfile << EOF
.tables
EOF

echo "REQUESTED '$table' TABLE."
echo "HEADER: "

sqlite3 $dbfile << EOF || exit 1
.headers on
.mode column
select *  from $table LIMIT 3;
EOF

echo "PICKING FIELDS: $fields"

sqlite3 $dbfile << EOF || exit 1
.mode csv
.separator  $sep
.headers on
.out $outputto
select $fields from $table;
EOF


echo "saved to $outputto"
