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
Reevaluates the quality of called SNPs (provided in VCF format). Returns:

- the reliability of a SNP (BOD score) based on [the described geometric method](http://dx.doi.org/10.1016%2Fj.ygeno.2012.12.001)

- per-nucleotide coverage of the reference around SNPs provided 
(+/- floor(1.5 * read length) around the SNP, tab-separated file, 
with blocks delimited by `\t|\t` ):

 - generally (100% and 90% identity)
 - for reads covering the SNP locus (100% and 90% identity)
 - location of centres of the reads covering the SNP locus


Sample files:
------------
The sample files for running the program are in the folder `testcase/`.


Additional information:
------------
A Python script (`src/plotCoverage.py`) is provided to produce the plots of the coverage.
An earlier R script (`src/rscript.R`) is outdated and is not compatible with the current version.

Contact:
--------
If you have any questions/concerns please drop an email to:
    fritz.sedlazeck@univie.ac.at
    d.lituiev@gmail.com
