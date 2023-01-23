# SsuisChara
SsuisChara is a tool for analyze the genome characteristics of Streptococcus suis using assemblied genome sequence for:
 * species verification -- Through align the 16S rRNA sequence, verify whether the input genomic sequence belongs to Streptococcus suis species or not.
 * Serotype prediction -- The capsular polysaccharide locus, of 29 classic streptococcus suis serotype were used as reference to find the best match in input genomic sequence to determine the serotype.
 * MLST sequence type -- Through screen 7 housekeeping genes of Streptococcus suis, determine the Sequence Type (ST) of input genomic sequence and allele number of every housekeeping gene
 * Virulence associated factors (vafs) screen -- Total 111 vafs of S. suis were collected from published papers and established as database to screen there presence and absence in input genomic sequence. 53 vafs distributed in accessory genome of S. suis, 58 vafs distributed in core genome of S. suis, two screen mode are provided, "concise" and "full", "concise" mode was set as default, only screen 53 vafs in accessory genome, "full" mode could screen all vafs. 
 * Clustermap of vafs generate -- A heatmap based on the presence and absence of vafs were generated and clustered to visualize the vafs prevalence.
 * Human infection potential -- Several vafs are found associated in previous study, there prevalence ratio in human S. suis isolated were used as weight of each vafs, we use the prevalence of these vafs to predict the human infection potential of each source isolate of input gneomic sequence.
 * antimicrobial resistance determinants -- Aminoglycoside, macrolide, and tetracycline are major class of antimicrobial drugs resist by S. suis, we screened the known amrgs resist these drugs to determine the antimicrobial resistance level of each source isolate of input gneomic sequence.

# External Dependencies
BLAST+
prodigal

# Usage
Put the python script and database folder into the folder contains your sequence file

```
SsuisChara [-i] [-o] [-r] [-t] [-v]
Input and Output:
  -i, --input             Input FASTA file
  -o, --output            Output file
Parameters:
  -t, --threads           Threads to use for BLAST searches
  --min_gene_cov          Minimum percentage coverage to consider a single gene complete. [default: 80.0%]
  --min_gene_id           Minimum percentage identity to consider a single gene complete. [default: 70.0%]
  --vf_screen_mode        The virulence factor screen mode, two modes, "concise" and "full" were provided. "concise" was set as default
  --no_heat_map           Suppress output of heatmap file
  -v, --version           Show version number and exit
```
# Quick usage
``` Python
python SsuisChara.py -i *.fasta 
```
