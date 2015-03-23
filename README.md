Quick start:
------------
    git clone --recursive https://github.com/Cibiv/BODscore

    cd BODscore
    bash ./install.sh

Functionality:
----------
Re-evaluates the quality of called SNPs (provided in VCF format). Returns:

- the reliability of a SNP (BOD score) based on [the described geometric method](http://dx.doi.org/10.1016%2Fj.ygeno.2012.12.001)

- per-nucleotide coverage of the reference around SNPs is returned either as:
  
  + a tab-delimited file
   (coverage for the range of `+/- floor(1.5 * read length)` around the SNP, 
   tab-separated file, with blocks delimited by `\t|\t` ):

  + an SQLite3 database with coverage written as blob columns (a faster option).
  The database can be read with the provided  Python `ReadCoverage.py` class 
  and related scripts

- following coverage data is returned within `QdrArray` objects:
  + accounting for perfectly matching reads with 100% identity (`HI`) 
    and loosely matching, between 90% and 100% identity (`LO`)
  + for forward and reverse strands

- three `QdrArray` objects with coverage profiles are returned:
  + for all reads in the range (`totCov` column in the SQLite database)
  + for reads covering only the locus proper (`snpCov` column)
  + location of centres of the reads covering the SNP locus (`alnCtr` column)

- the coverage data is written in a `sampleX__coverage` table with a composite key of 
`(contig, pos)`, where `contig` is of `TEXT` type and `pos` (position) of `INT` type.

Test sample files:
------------
Test sample files for running the program are in the folder `testcase/`.


Additional information:
------------
A Python script (`src/plot_coverage_sqlite.py`) is provided to produce the plots of the coverage.
An earlier R script (`src/rscript.R`) is outdated and is not compatible with the current version.


Contact:
--------
If you have any questions/concerns please drop an email to:

    fritz.sedlazeck@gmail.com

    d.lituiev@gmail.com
