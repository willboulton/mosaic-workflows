#!/bin/bash -e

DIR=/gpfs/data/mock_lab/will/mosaic_tmp/binning_dir
SAMPLE=$1

DEPTH=$DIR/$SAMPLE/$SAMPLE.a.depth.txt
DEPTH2=$DIR/$SAMPLE/$SAMPLE.a.depth2.txt

if [ ! -f "$DEPTH2" ]; then
    mv $DEPTH $DEPTH2
fi

printf "ID\tAvgFold\n" > $DEPTH

tail -n +2 $DEPTH2 | sed 's|_|\t|2' | cut -f 1,3 | uniq >> $DEPTH