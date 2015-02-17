#!/bin/bash

cp libs/sqlite3pp_src_CMakeLists.txt  libs/sqlite3pp/src/CMakeLists.txt

unzip "libs/sqlite-amalgamation-3080802.zip" -d "libs/sqlite3pp/src/" || exit 1;

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
