vshape="../bin/vshape-1.0.0/vshape"
vcf="test.sort.vcf"
bam="short.bam"
REFERENCE="Ref.fa"
FIRST_N_ENTRIES=5 # only this number of entries will be read. normally just omit this flag or set to 0

DATABASE="test.db"

$vshape -v -b $FIRST_N_ENTRIES -x mitochondria -x chloroplast -l 50 -r $REFERENCE  -d  $DATABASE -s $vcf -q $bam -t "testsample" || exit 1

echo -e "output has been saved to:\t$DATABASE" >&2
read -r -d '' info <<EOF
now you can :
(1) visualize the local coverage around the SNPs with python scripts, e.g.:
\tvisualise.sh
(2) dump the content to a csv / tab file. See:
\tdbtocsv.sh
(3) use database viewers (e.g. http://sqlitebrowser.org/) to see the database content
EOF

echo -e "$info" >&2
