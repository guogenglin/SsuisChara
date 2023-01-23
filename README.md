# SsuisChara
SsuisChara is a tool for analyze the genome characteristics of Streptococcus suis using assemblied genome sequence for:
 * species verification -- Through align the 16S rRNA sequence, verify whether the input genomic sequence belongs to Streptococcus suis species or not
 * Serotype prediction -- The capsular polysaccharide locus, of 29 classic streptococcus suis serotype were used as reference to find the best match in input genomic sequence to determine the serotype.
 * MLST sequence type -- Through screen 7 housekeeping genes of Streptococcus suis, determine the Sequence Type (ST) of input genomic sequence and allele number of every housekeeping gene
 * Virulence associated factors (vafs) screen -- Total 111 vafs of S. suis were collected from published papers and established as database to screen there presence and absence in input genomic sequence
 * Clustermap of vafs generate -- A heatmap based on the presence and absence of vafs were generated and clustered to visualize the vafs prevalence
 * Human infection potential -- Several vafs are found associated in previous study, there prevalence ratio in human S. suis isolated were used as weight of each vafs, we use the prevalence of these vafs to predict the human infection potential of each source isolate of input gneomic sequence
 * antimicrobial resistance determinants -- Aminoglycoside, macrolide, and tetracycline are major class of antimicrobial drugs resist by S. suis, we screened the known amrgs resist these drugs to determine the antimicrobial resistance level of each source isolate of input gneomic sequence
