#!/bin/bash

for f in `ls $1`; do 
    for g in `ls $1`; do 
        printf "${f##*/}\t${g##*/}\t"
        comm -12 \
            <(grep '>' $1/$f | tr -d '>' |sort) \
            <(grep '>' $1/$g | tr -d '>' | sort) | wc -l
    done 
done > $2