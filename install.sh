#!/bin/bash

cp libs/sqlite3pp_src_CMakeLists.txt  libs/sqlite3pp/src/CMakeLists.txt

FLAG_SQLITE_INCLUDED=1;

if [ ! $FLAG_SQLITE_INCLUDED  ]
then 
    echo "Downloading SQLite3 base..."
    wget  -N -q -A.zip "http://www.sqlite.org/2015/sqlite-amalgamation-3080802.zip" -P "./libs" || exit 1;
    cp "libs/sqlite3pp_src_sqlite_CMakeLists.txt" "libs/sqlite3pp/src/sqlite-amalgamation-3080802/CMakeLists.txt"
    unzip "libs/sqlite-amalgamation-3080802.zip" -d "libs/sqlite3pp/src/" || exit 1;
else
    echo "====================================================="
fi

mkdir build
cd build || exit 1;
rm ./* -rf

cmake .. || exit 1;

make || exit 1;

cd ../bin || ( echo "hmm. something went wrong... \n No 'bin' directory has been produced."; exit 1;)

ls -l

echo "======================================================="
echo "go to the directory './vshape-#.#.#' and run './vshape'"
echo " good luck! "
