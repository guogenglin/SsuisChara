#!/bin/bash

# profiles
curl https://rest.pubmlst.org/db/pubmlst_ssuis_seqdef/schemes/1/profiles_csv \
 -o profiles.tsv
 
# loci
for locus in aroA cpn60 dpr gki mutS recA thrA
do
    curl "https://rest.pubmlst.org/db/pubmlst_ssuis_seqdef/loci/${locus}/alleles_fasta" \
     -o ${locus}.fasta
done

# convert fasta to one line
for f in *.fasta; do
    seqtk seq "$f" > temp.fasta && mv temp.fasta "$f"
done