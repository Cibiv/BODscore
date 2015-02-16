Quick start:
------------

    mkdir -p build/release
    cd build/release/
    cmake ../..
    make
    cd ../../bin/vshape-*
    ./vshape


Functionality:
----------
Re-evaluates the quality of called SNPs (provided in VCF format). Returns:

- the reliability of a SNP (BOD score) based on [the described geometric method](http://dx.doi.org/10.1016%2Fj.ygeno.2012.12.001)

- per-nucleotide coverage of the reference around SNPs is returned either as:
  
  + a tab-delimited file
   (coverage for the range of +/- floor(1.5 * read length) around the SNP, 
   tab-separated file, with blocks delimited by `\t|\t` ):

  + an SQLite3 database with coverage written as blob columns (a faster option).
  The database can be read with the provided  Python `ReadCoverage.py` class 
  and related scripts

- following coverage data is returned within one block:
  + accounting for perfectly matching reads with 100% identity ('`HI`') 
    and loosely matching, between 90% and 100% identity ('`LO`')
  + for forward and reverse strands

- three blocks of coverage profiles are returned:
  + for all reads in the range (`totCov` column in the SQLite database)
  + for reads covering only the locus proper (`snpCov` column)
  + location of centres of the reads covering the SNP locus (`alnCtr` column)

- each chromosome is written in a separate table within the database:

   `sample123__coverage_1` -- chr1,

   `sample123__coverage_2` -- chr2 etc.


Test sample files:
------------
Test sample files for running the program are in the folder `testcase/`.


Additional information:
------------
A Python script (`src/plotCoverage.py`) is provided to produce the plots of the coverage.
An earlier R script (`src/rscript.R`) is outdated and is not compatible with the current version.


Contact:
--------
If you have any questions/concerns please drop an email to:

    fritz.sedlazeck@gmail.com

    d.lituiev@gmail.com
