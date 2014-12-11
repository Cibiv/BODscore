Quick start:
------------

mkdir -p build/release
cd build/release/
cmake ../..
make
cd ../../bin/vshape-*
./vshape


Options:
----------
There are 2 options implemented:

1: Summarize identity
	Computes the identity distribution along mapped reads.
		
2: Reevaluates the called SNPs
	Uses the described geometric method to compute the reliability of a SNP



Additional information:
------------
A R script is provided (src/rscript.R) to produce the plots of the coverage


Contact:
------------
If you have any questions/concerns please drop me an email:
fritz.sedlazeck@univie.ac.at
