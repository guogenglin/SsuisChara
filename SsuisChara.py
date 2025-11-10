# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 10:41:14 2025

@author: Genglin Guo
"""

import argparse
import sys
import pathlib
import multiprocessing
import subprocess
import time
from Bio import SeqIO
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt

__version__ = '2.0'

def get_argument():
    # Parsers
    parser = argparse.ArgumentParser(description = 'SsuisChara', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_group_1 = parser.add_argument_group('Input and Output')
    parser_group_2 = parser.add_argument_group('Parameters')

    # Input and output
    parser_group_1.add_argument('-i', '--input', required = True, nargs = '+', type = str, 
                                help = 'Input FASTA file')
    parser_group_1.add_argument('-o', '--output', required = False, type = str, default = 'SsuisChara_output.txt',
                              help = 'Output file')

    # Parameters
    parser_group_2.add_argument('-t', '--threads', required = False, type = int, default = min(multiprocessing.cpu_count(), 4), 
                        help = 'Threads to use for BLAST searches')
    parser_group_2.add_argument('--min_gene_cov', required = False, type = float, default = 80.0, 
                               help = 'Minimum percentage coverage to consider a single gene complete. [default: 80.0%]')
    parser_group_2.add_argument('--min_gene_id', required = False, type = float, default = 70.0, 
                               help = 'Minimum percentage identity to consider a single gene complete. [default: 70.0%]')
    parser_group_2.add_argument('--vf_screen_mode', required = False, type = str, default = 'concise', 
                                help = 'The virulence factor screen mode, two modes, "concise" and "full" were provided. "concise" was set as default.')
    parser_group_2.add_argument('--vf_location_details', action='store_true', 
                                help = 'Generate a table contain the location of every screened virulence factors in the input sequence.')
    parser_group_2.add_argument('--heatmap', action='store_true', help = 'Generate a heatmap file to visualize the prevalence of VFs')
    parser_group_2.add_argument('-v', '--version', action = 'version', version = 'SsuisChara v' + __version__, 
                        help = 'Show version number and exit')
    return parser

def check_dependencies():
    # Checks dependencies are available
    dependencies = ['makeblastdb', 'blastn', 'blastx', 'prodigal']
    for i in dependencies:
        try:
            subprocess.check_call(['which', i], stdout = subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print('Error: could not find %s tool' % i, file=sys.stderr)
            sys.exit(1)

def get_best_species_result(inputfile, refpath, threads):
    # Perform blast, find and align the 16S sequence with the inputfile, pending the species with a threshold 97%
    repa = pathlib.Path(refpath).resolve() / 'Ssuis16s.fas'
    inpa = pathlib.Path(inputfile).resolve()
    blast_hits = run_blast(inpa, repa, threads, 'n')
    species = ''
    for hit in blast_hits:
        if hit.length < 1000:
            continue
        if hit.pident >= 97.0:
            species = 'Streptococcus suis'
    if not species:
        species = 'NA'
    return species

def run_blast(inpa, repa, threads, setting):
    # Do blast, iterator the result to a list
    blast_hits = []
    if setting == 'mlst':
        command = ['blastn', '-query', inpa, '-subject', repa, '-num_threads', str(threads), '-evalue', '0.00001', '-perc_identity', '99.5', 
                   '-outfmt', '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    elif setting == 'vf':
        command = ['tblastn', '-query', inpa, '-subject', repa, '-num_threads', str(threads), '-outfmt', 
                   '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    else:
        command = ['blastn', '-query', inpa, '-subject', repa, '-num_threads', str(threads), '-outfmt', 
                   '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    
    process = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    try:
        for line in out.strip().split('\n'):
            blast_hits.append(BlastResult(line))
    except:
        blast_hits = []
    return blast_hits

class BlastResult(object):
    # Handle the BLAST output
    def __init__(self, hit_string):
        parts = hit_string.split('\t')
        self.qseqid = parts[0]
        self.sseqid = parts[1]
        self.qstart = int(parts[2]) - 1
        self.qend = int(parts[3])
        self.sstart = int(parts[4]) - 1
        self.send = int(parts[5])
        if self.sstart <= self.send:
            self.strand = '+'
        else:
            self.sstart, self.send = self.send, self.sstart
            self.strand = '-'
        self.length = int(parts[8])
        self.pident = float(parts[9])
        self.query_cov = 100.0 * len(parts[11]) / float(parts[10])

def get_serotype(inputfile, refpath, threads):
    # Find the best serotype of inputfile by sequence alignment with 33 reference cps locus sequences
    repa = pathlib.Path(refpath).resolve() / 'Streptococcus_suis_cps_locus.fasta'
    inpa = pathlib.Path(inputfile).resolve()
    blast_hits = run_blast(repa, inpa, threads, 'n')
    refs, blast_result = {}, {}
    for serotype in SeqIO.parse(repa, 'fasta'):
        refs[serotype.name] = serotype.seq
    for hit in blast_hits:
        if not hit.qseqid in blast_result:
            blast_result[hit.qseqid] = []
        if hit.qend > hit.qstart:
            blast_result[hit.qseqid].append((hit.qstart, hit.qend, hit.pident))
        else:
            blast_result[hit.qseqid].append((hit.qend, hit.qstart, hit.pident))
    simplified_result = simplify(blast_result)
    best_serotype = ''
    best_coverage = 0.0
    best_identity = 0.0
    for key, value in simplified_result.items():
        length = get_total_length(value)
        length_identity = get_total_length_identity(value)
        coverage = 100.0 * length / len(refs[key])
        identity = 100.0 * length_identity / len(refs[key])
        if coverage > best_coverage:
            best_serotype = key
            best_coverage = coverage
            best_identity = identity
        elif coverage == best_coverage and identity > best_identity:
            best_serotype = key
            best_identity = identity
    if not best_serotype:
        best_serotype = 'NA'
    if best_coverage <= 95.0:
        best_serotype += '?'
    if best_identity > 100.00:
        best_identity = 100.00
    return best_serotype, round(best_coverage, 2), round(best_identity, 2)

def simplify(blast_result):
    # Blast result may have multiple result and overlap, simplify the result, and convert the percent_identity to length_identity
    simplified_result = {}
    for key, value in blast_result.items():
        simplified_list = []
        ranks = sorted(value, key = lambda v: v[0])
        start = 0
        end = 0
        length_ident = 0
        for rank in ranks:
            if rank[0] <= end:
                if rank[1] > end:
                    end = rank[1]
                    length_ident += (rank[1] - rank[0]) * rank[2] / 100
                else:
                    continue
            else:
                if end == 0:
                    length_ident += (rank[1] - rank[0]) * rank[2] / 100
                    start, end = rank[0], rank[1]
                else:
                    length_ident += (rank[1] - rank[0]) * rank[2] / 100
                    simplified_list.append((start, end, length_ident))
                    start, end, length_ident = rank[0], rank[1], 0
        simplified_list.append((start, end, length_ident))
        simplified_result[key] = simplified_list
    return simplified_result

def get_total_length(value):
    # Sum the length of all non-redundant blast hit of one single cps locus sequence
    cover_length = 0
    for hit in value:
        length = hit[1] - hit[0]
        cover_length += length
    return cover_length

def get_total_length_identity(value):
    # Sum the identified length of all non-redundant blast hit of one single cps locus sequence
    length_identity = 0
    for hit in value:
        length_identity += hit[2]
    return length_identity

def get_best_mlst_result(inputfile, refpath, mlst_database, labels, threads):
    # Perform blast, find all allele numbers and compare with reference, collect best ST match
    inpa = pathlib.Path(inputfile).resolve()
    best_match = []
    best_match_in_ref = []
    best_ST = ''
    for i in labels:
        filename = i + '.fas'
        repath = refpath / 'MLST' / filename
        repa = pathlib.Path(repath).resolve()
        blast_hits = run_blast(repa, inpa, threads, 'mlst')
        best_allele = ''
        best_cov = 0.0
        best_pident = 0.0
        for hit in blast_hits:
            if hit.length < 100:
                continue
            if hit.pident >= best_pident and hit.query_cov >= best_cov:
                best_pident = hit.pident
                best_cov = hit.query_cov
                best_allele = hit.qseqid.split('_')[-1]
        if best_cov == 100.0 and best_pident == 100:
            best_match.append(best_allele)
            best_match_in_ref.append(best_allele)
        else:
            best_match_in_ref.append(best_allele)
            best_allele += '?'
            best_match.append(best_allele)
    for key, value in mlst_database.items():
        if value == best_match_in_ref:
            best_ST = key
    if any('?' in item for item in best_match) and best_ST:
        best_ST += '*'
    if not best_ST:
        best_ST = 'NA'
    return best_ST, best_match

def process_mlst_reference(refpath):
    # Generate a dict for mlst, ST : alleles and collect the mlst labels
    database = {}
    labels = []
    profilepath = pathlib.Path(refpath).resolve() / 'MLST' / 'MLST_profiles'
    with open(profilepath, 'rt') as file:
        for line in file:
            profile = line.strip('\n').split('\t')
            if not labels:
                for i in profile[1:]:
                    labels.append(i)
            else:
                ST = profile[0]
                alleles = profile[1:]
                database[ST] = alleles
    return database, labels

def get_amrg_result(inputfile, refpath, threads, min_gene_id, min_gene_cov):
    # Find the aminoglycoside, macrolide, and tetracycline AMRGs in the input genomic sequence
    AMRG_level = 0
    amrg_dict = {}
    inpa = pathlib.Path(inputfile).resolve()
    for i in ['aminoglycoside', 'macrolide', 'tetracycline']:
        file_name = i + '.fas'
        repath = refpath / 'AMR' / file_name
        repa = pathlib.Path(repath).resolve()
        all_hit_result = []
        blast_hits = run_blast(repa, inpa, threads, 'n')
        if blast_hits:
            for hit in blast_hits:
                if hit.length < 100:
                    continue
                if hit.pident >= min_gene_id and hit.query_cov >= min_gene_cov:
                    best_match = hit.qseqid.split('_')[0]
                    if best_match not in all_hit_result:
                        all_hit_result.append(best_match)
        if all_hit_result:
            AMRG_level += 1
        else:
            all_hit_result.append('NA')
        amrg_dict[i] = ';'.join(all_hit_result)
    return AMRG_level, amrg_dict['aminoglycoside'], amrg_dict['macrolide'], amrg_dict['tetracycline']

def get_vfs_result(inputfile, refpath, threads, mode, min_gene_id, min_gene_cov, location_details):
    # Screen the VF genes in input genomic sequence
    vf_genes = []
    vf_details = dict()
    repath = refpath / 'VFs'
    inpa = pathlib.Path(inputfile).resolve()
    if mode == 'concise':
        repa = pathlib.Path(repath).resolve() / 'vfs_in_accessory_genome.fasta'
    elif mode == 'full':
        repa = redirection(repath)
    blast_hits = run_blast(repa, inpa, threads, 'vf')
    if blast_hits:
        for hit in blast_hits:
            best_match = ''
            if hit.length < 100:
                continue
            if hit.pident >= min_gene_id and hit.query_cov >= min_gene_cov:
                best_match = hit.qseqid
            if best_match and best_match not in vf_genes:
                vf_genes.append(best_match)
                if best_match not in vf_details:
                    gene_code = hit.sseqid + '_' + str(hit.sstart) + '_' + str(hit.send) + '_' + hit.strand
                    vf_details[best_match] = gene_code
                else:
                    vf_details[best_match] += '; '
                    vf_details[best_match] += gene_code
    vfs_counts = len(vf_genes)
    headers = find_all_vfs_name(repa)
    if location_details:
        generate_vf_matrix_output(inputfile, vf_details, headers)
    if mode == 'full':
        pathlib.Path(repa).unlink()
    return vfs_counts, vf_genes

def redirection(repath):
    # Based on the given command, redirect the vf gene database if necessary
    newpath = pathlib.Path(repath).resolve() / 'vfs_in_pan_genome.fasta'
    corepa = pathlib.Path(repath).resolve() / 'vfs_in_core_genome.fasta'
    accepa = pathlib.Path(repath).resolve() / 'vfs_in_accessory_genome.fasta'
    core = open(corepa, 'rt')
    accessory = open(accepa, 'rt')
    with open(newpath, 'at') as file:
        file.write(accessory.read())
        file.write(core.read())
    core.close()
    accessory.close()
    return newpath

def find_all_vfs_name(repa):
    # Generate a blank output table file
    headers = ['Isolate']
    for gene in SeqIO.parse(repa, 'fasta'):
        headers.append(gene.name)
    return headers

def generate_vf_matrix_output(inputfile, vf_details, headers):
    # Generate output
    if not pathlib.Path('vf_matrix.txt').is_file():
        with open('vf_matrix.txt', 'wt') as file:
            file.write('\t'.join(headers))
            file.write('\n')
    line = [inputfile]
    for vf in headers[1:]:
        if vf in vf_details:
            line.append(vf_details[vf])
        else:
            line.append('NA')
    with open('vf_matrix.txt', 'at') as file:
        file.write('\t'.join(line))
        file.write('\n')  

def calculate_zoonotic_potential(vf_genes):
    # Based on the vf genes' presence and absence, determine the zoonotic potential
    score_for_vfs = {
        'cbp40omp40' : 3.003,
        'cps2B' : 1.877,
        'cps2E' : 2.402,
        'cps2F' : 2.402,
        'cps2G' : 2.402,
        'cps2J' : 2.477,
        'cps2L' : 2.102,
        'Fhb_1' : 10.210,
        'Fhb_2' : 3.754,
        'Hhly3' : 7.132,
        'hylA' : 3.078,
        'ldeS' : 3.979,
        'lgdE' : 4.955,
        'MRP' : 5.330,
        'NeuB' : 2.027,
        'NisK' : 9.234,
        'NisR' : 9.309,
        'PnuC' : 4.505,
        'refA' : 2.327,
        'RggTR' : 4.354,
        'sbp1' : 2.853,
        'sbp2' : 3.003,
        'sly' : 1.8018,
        'SP1' : 1.577,
        'Tran' : 2.027,
        'Zmp' : 1.877
        }
    zoonotic_score = 0.0
    human_infection_potential = ''
    for gene in vf_genes:
        if gene in score_for_vfs:
            zoonotic_score += score_for_vfs[gene]
        else:
            continue
    zoonotic_score = round(zoonotic_score, 3)
    if zoonotic_score >= 70.0:
        human_infection_potential = 'high'
    elif 30.0 <= zoonotic_score < 70.0:
        human_infection_potential = 'medium'
    else:
        human_infection_potential = 'low'
    return human_infection_potential, zoonotic_score

def generate_output(output, labels, refpath, mode):
    # Generate a blank output table file
    repath = refpath / 'VFs'
    if mode == 'concise':
        repa = pathlib.Path(repath).resolve() / 'vfs_in_accessory_genome.fasta'
    elif mode == 'full':
        repa = redirection(repath)
    headers = find_all_vfs_name(repa)
    total_vfs = str(len(headers) - 1)
    if pathlib.Path(output).is_file():
        return
    vfs = 'vfs counts (n/' + total_vfs + ')'
    headers = ['Isolate', 'Species', 'Serotype', 'Coverage', 'Identity', 'ST', vfs, 'human infection potential', 
               'zoonotic_score', 'AMRG_level', 'aminoglycoside', 'macrolide', 'tetracycline']
    for i in labels:
        headers.append(i)
    with open(output, 'wt') as file:
        file.write('\t'.join(headers))
        file.write('\n')

def positive_output(output, inputfile, species, serotype, sero_coverage, sero_identity, best_ST, VFs_counts, vf_genes, human_infection_potential, 
       zoonotic_score, AMRG_level, amino, macro, tetra, best_mlst_match):
    # Generate output
    simple_output = inputfile + ' : ' + species + ' ' + 'serotype' + ' ' + serotype + ' ' + 'ST' + ' ' + best_ST
    line = [inputfile, species, serotype, str(sero_coverage), str(sero_identity), best_ST, str(VFs_counts), human_infection_potential, str(zoonotic_score), 
            str(AMRG_level), amino, macro, tetra]
    for i in best_mlst_match:
        line.append(i)
    print(simple_output)
    with open(output, 'at') as file:
        file.write('\t'.join(line))
        file.write('\n') 

def negative_output(output, inputfile):
    # Generate output
    simple_output = inputfile + ' : ' + 'No Streptococcus suis'
    line = [inputfile, 'NA']
    print(simple_output)
    with open(output, 'at') as file:
        file.write('\t'.join(line))
        file.write('\n') 

def generate_heatmap(all_input_names, refpath, all_vfs, mode):
    # Generate a virulence gene presence and absence heatmap
    repath = refpath / 'VFs'
    if mode == 'concise':
        repa = pathlib.Path(repath).resolve() / 'vfs_in_accessory_genome.fasta'
    elif mode == 'full':
        repa = redirection(repath)
    vf_names = []
    with open(repa, 'rt') as file:
        for line in file:
            if line.startswith('>'):
                line = line.strip('\n')
                vf_names.append(line[1:])
            else:
                continue
    vf_matrix = []
    for all_vf in all_vfs:
        vf_line = []
        for vf_name in vf_names:
            if vf_name in all_vf:
                vf_line.append(1)
            else:
                vf_line.append(0)
        vf_matrix.append(vf_line)
    matrix = pd.DataFrame(vf_matrix)
    matrix.columns = vf_names
    matrix.index = all_input_names
    figx = len(vf_names) / 3
    figy = len(all_input_names)
    sns.clustermap(matrix, cmap = 'Set2_r', linewidth=1, col_cluster = False, figsize = (figx, figy))
    savepath = pathlib.Path(refpath).parent
    savepath = savepath / 'heatmap.png'
    plt.savefig(savepath, dpi = 600)
    if mode == 'full':
        pathlib.Path(repa).unlink()

def main():
    print('If you have any questions or suggestions for SsuisChara, please contact Genglin Guo, e-mail: 2019207025@njau.edu.cn')
    starttime = time.perf_counter()
    # Initialize
    args = get_argument().parse_args()
    check_dependencies()
    # Prepare for potential heatmap generate
    if args.heatmap:
        all_input_names = []
        all_vfs = []
    refpath = pathlib.Path(__file__).resolve().parent / 'database'
    # process the mlst database
    mlst_database, labels = process_mlst_reference(refpath)
    # Run this pipeline for each single input genome
    for inputfile in args.input:
        if args.heatmap:
            all_input_names.append(inputfile)
        species = get_best_species_result(inputfile, refpath, args.threads)
        if species == 'Streptococcus suis':
            serotype, sero_coverage, sero_identity = get_serotype(inputfile, refpath, args.threads)
            best_ST, best_mlst_match = get_best_mlst_result(inputfile, refpath, mlst_database, labels, args.threads)
            AMRG_level, aminoglycoside, macrolide, tetracycline = get_amrg_result(inputfile, refpath, args.threads, args.min_gene_id, args.min_gene_cov)
            vfs_counts, vf_genes = get_vfs_result(inputfile, refpath, args.threads, args.vf_screen_mode, args.min_gene_id, args.min_gene_cov, args.vf_location_details)
            human_infection_potential, zoonotic_score = calculate_zoonotic_potential(vf_genes)
            if args.heatmap:
                all_vfs.append(vf_genes)
    # Generate output
            generate_output(args.output, labels, refpath, args.vf_screen_mode)
            positive_output(args.output, inputfile, species, serotype, sero_coverage, sero_identity, best_ST, vfs_counts, vf_genes, human_infection_potential, 
               zoonotic_score, AMRG_level, aminoglycoside, macrolide, tetracycline, best_mlst_match)
        else:
            generate_output(args.output, labels, refpath, args.vf_screen_mode)
            negative_output(args.output, inputfile)
    # Generate heatmap
    if  args.heatmap:
        generate_heatmap(all_input_names, refpath, all_vfs, args.vf_screen_mode)
    # Total time count
    endtime = time.perf_counter() - starttime
    per_genome_time = endtime / len(args.input)
    print('{:.1f}h{:.1f}m{:.1f}s for one genome'.format(per_genome_time // 3600, per_genome_time % 3600 // 60, per_genome_time % 60))
    print('Total time consumed : {:.1f}h{:.1f}m{:.1f}s'.format(endtime // 3600, endtime % 3600 // 60, endtime % 60))
   
if __name__ == '__main__':

    main()

