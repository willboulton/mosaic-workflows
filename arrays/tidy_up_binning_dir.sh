#!/bin/bash

SAMPLE=$1
BD=/gpfs/data/mock_lab/will/mosaic_tmp/binning_dir

echo "tidying $SAMPLE"
cd $BD/$SAMPLE


echo "removing files from zip file"
find $BD/$SAMPLE/$SAMPLE.a.fna -delete
find $BD/$SAMPLE/$SAMPLE.a.faa -delete
find $BD/$SAMPLE/$SAMPLE.a.map.txt -delete
find $BD/$SAMPLE/$SAMPLE.a.depth.txt -delete
find $BD/$SAMPLE/$SAMPLE.a.pfam.txt -delete
find $BD/$SAMPLE/$SAMPLE.*.gff -delete
find $BD/$SAMPLE/$SAMPLE.bins.faa.* -delete
find $BD/$SAMPLE/$SAMPLE.bins.faa -delete
find $BD/$SAMPLE/$SAMPLE.bins.fna -delete

echo "removing semibin data"
find $BD/$SAMPLE/semibin/data.csv -delete
find $BD/$SAMPLE/semibin/data_split.csv -delete

echo "removing vamb data"
find $BD/$SAMPLE/vamb/composition.npz -delete  
find $BD/$SAMPLE/vamb/contignames -delete
find $BD/$SAMPLE/vamb/latent.npz -delete
find $BD/$SAMPLE/vamb/lengths.npz -delete
find $BD/$SAMPLE/vamb/mask.npz -delete
find $BD/$SAMPLE/vamb/model.pt -delete

echo "removing directories"
rm -rf $BD/$SAMPLE/bins_checkm/protein_files
rm -rf $BD/$SAMPLE/bins_checkm/diamond_output
rm -rf $BD/$SAMPLE/semibin/output_bins
rm -rf $BD/$SAMPLE/vamb/bins
rm -rf $BD/$SAMPLE/bins_das/das_DASTool_bins
rm -rf $BD/$SAMPLE/ebins_eukcc/refine_workdir

echo "make metabat tsv"
grep -H ">" $BD/$SAMPLE/metabat/* | \
    sed -e 's|.fa:>|\t|g' -e 's|metabat/||g' > $BD/$SAMPLE/metabat/metabat_c2b.tsv && \
    rm -f $BD/$SAMPLE/metabat/*.fa

echo "zipping files"
tar czf $BD/$SAMPLE/bins.tar.gz $BD/$SAMPLE/bins && rm -rf $BD/$SAMPLE/bins
tar czf $BD/$SAMPLE/ebins.tar.gz $BD/$SAMPLE/ebins && rm -rf $BD/$SAMPLE/ebins

echo "done"