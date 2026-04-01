#!/bin/bash

# three class
for amr in aminoglycoside macrolide tetracycline
do
    curl "https://bitbucket.org/genomicepidemiology/resfinder_db/raw/master/${amr}.fsa" \
     -o ${amr}.fas
done

# convert fasta to one line
for f in *.fas; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"
done