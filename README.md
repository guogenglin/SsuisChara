# SsuisChara
SsuisChara is a tool for integrally analyze the genome characteristics of *Streptococcus suis* using assemblied genome sequence for:
 * Species verification
 * Serotype prediction
 * MLST sequence type
 * Virulence associated factors (vafs) screen
 * Human infection potential
 * Antimicrobial resistance determinants

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
# Species verification
  The 16S rRNA sequence of *Streptococcus suis* were obtained from public database, through align the 16S rRNA sequence with input sequence, the 16S sequence in input sequence could be found. If the identity of 16S rRNA sequence lower than 97%, we consider the input sequence do not belong to *Streptococcus suis* species, otherwise, we consider the input sequence belong to *Streptococcus suis* species.
# Serotype prediction
  Initially, serotyping was based on serological tests, subsequently, the serotyping determine locus of many bacteria were found in the genome, and a lot of molecular serotyping method were developed based on the difference of serotyping determine locus, such as multiplex PCR. Now the number of acquired sequenced bacteria genome in public databased are fastly increasing, allow us to explore a full locus alignment method to high throughput and precisely determine serotype *in silico*.
  
  For *Streptococcus suis*, the serotype determine locus is capsular polysaccharide (cps) locus. In this analysis, the cps locus of 29 classic *Streptococcus suis* serotype were collected and used as reference to find the best match in input genome to determine the serotype. Dispite a lot of novel capsular polysaccharide loci (NCL) were found, however, the importance and prevalence of these serotype as limited, so we haven't include them, however, we may upload in the future.
  
  To aviod the uncorrected prediction, we displayed the coverage and identity of predicted serotype, if it lower than the threshold (we set as 95% prelimitarily), a "?" will be added to the end of output string of predicted serotype.
# MLST sequence type
  The MLST database were obtained from [pubMLST](https://pubmlst.org/), update: 2023-01-12
  
  Seven housekeeping genes, *aroA*, *cpn60*, *dpr*, *gki*, *mutS*, *recA*, and *thrA*, of *Streptococcus suis*, were screened in input genome and return there allele number, or closest allele number, then determine there Sequence Type (ST), if not every allele number are exact match, a "?" will be added to the end of output string of predicted ST.
# Virulence associated factors (vafs) screen
  Total 111 vafs of *S. suis* were collected from published papers and established as database to screen there presence and absence in input genome. 53 vafs distributed in accessory genome of S. suis, 58 vafs distributed in core genome of *S. suis*, two screen mode are provided, "concise" and "full", "concise" mode was set as default, only screen 53 vafs in accessory genome, "full" mode could screen all vafs.
  
  A heatmap based on the presence and absence of vafs were generated and clustered to visualize the vafs prevalence. If the user do not want to generate a heatmap, could use '''--no_heat_map''' command.
# Human infection potential
  Several vafs are found associated with human *S. suis* infection in previous study, there prevalence ratio in human *S. suis* isolated were used as weight of each vafs, we use the prevalence of these vafs to predict the human infection potential of each source isolate of input gneomic sequence.
  
  The summation of weights were named as "zoonotic score", if zoonotic_score >= 70.0, we give human_infection_potential as "high", if 30.0 <= zoonotic_score < 70.0, we give human_infection_potential as "medium", if zoonotic_score < 30.0, we give human_infection_potential as "low".
# Antimicrobial resistance determinants
  Aminoglycoside, macrolide, and tetracycline are major class of antimicrobial drugs resist by *S. suis*, we screened the known AMRGs resist these drugs to determine the antimicrobial resistance level of each source isolate of input gneomic sequence, the AMRG_level = The number of these 3 AMR drugs class covered by the AMRG screened in input genome. Then the screened AMRGs will be output in output file.
# Output
  The simlified result will generated in terminal, a detailed table were generated in work folder, also a heatmap.
