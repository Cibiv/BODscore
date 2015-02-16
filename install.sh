#!/bin/bash

cp libs/sqlite3pp_src_CMakeLists.txt  libs/sqlite3pp/src/CMakeLists.txt

mkdir build
cd build || exit 1;
rm ./* -rf

cmake ..

make

cd ../bin

ls -l

echo "======================================================="
echo "go to the directory './vshape-#.#.#' and run './vshape'"
echo " good luck! "
