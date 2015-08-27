vshape="../bin/vshape-1.0.0/vshape"
vcf="test.sort.vcf"
bam="short.bam"
REFERENCE="Ref.fa"
FIRST_N_ENTRIES=5 # only this number of entries will be read. normally just omit this flag or set to 0

TXTOUT="test.out.tab"
$vshape -v -b $FIRST_N_ENTRIES -x mitochondria -x chloroplast -l 50 -r $REFERENCE  -p  $TXTOUT -s $vcf -q $bam -t "testsample"
echo "output written to:\t$TXTOUT" >&2


TXTOUT_MINIMAL="test.out.min.tab"
echo -e "contig\tposition\tscore" > "$TXTOUT_MINIMAL"
cut -f1-3  test.out.tab > "$TXTOUT_MINIMAL"
echo "minimal output written to:\t$TXTOUT_MINIMAL" >&2
