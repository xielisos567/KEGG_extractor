#!/bin/python3

import os
import re
import sys
import glob
import argparse
import pandas as pd
from collections import defaultdict

def getSeq(fafile, genelist):
    '''
    给定fa文件，提取目标序列
    '''
    rank = 0
    out_seq = ''
    fa_file = open(fafile)
    for line in fa_file:
        if line.startswith('>') and rank==0:
            seq_id = line.strip()
            if seq_id.split(' ')[0][1:] == genelist or (re.search('protein_id=.*?\]',seq_id)!=None and re.search('protein_id=(.*?)\]',seq_id).group(1) == genelist):
                rank = 1
        elif line.startswith('>') and rank==1:
            break
        elif rank == 1:
            out_seq += line.strip()
        else:
            continue
    return seq_id, out_seq

def speciesExtract(species, species_info):
    '''
    提取物种信息
    '''
    print (species)
    try:
        try:
            species_file = pd.read_table(species, usecols=[0,7,8],header=None,engine='python')
        except:
            species_file = pd.read_table(species, encoding='unicode_escape', usecols=[0,7,8],header=None,engine='python')
    except:
        return species_info
    for i in range(species_file.shape[0]):
        #print (species_file.iloc[i])
        if str(species_file.iloc[i,2]) == 'nan':
            species_info[species_file.iloc[i,0]] = species_file.iloc[i,1]
        else:
            species_info[species_file.iloc[i,0]] = ' '.join([species_file.iloc[i,1], str(species_file.iloc[i,2])])
    return species_info

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='extract seq')
    parser.add_argument("-i", '--refdir', required=True, help='reference dir')
    parser.add_argument("-f", '--result', required=True, help="analysis result")
    parser.add_argument("-s", '--speciesdir', required=True, help='species info dir')
    parser.add_argument("-o", '--outdir', required=True, help='output dir')
    args = parser.parse_args()
    species = {}
    for file in glob.glob(f'{args.speciesdir}/*txt'):
        #print (file)
        species = speciesExtract(file, species)
    analysis_out = open(args.result)
    target_seq = defaultdict(list)
    ko_info = defaultdict(list)
    for line in analysis_out:
        line_info = re.split('\s{1,}',line.strip())
        genome_file = ''
        ref_info = line_info[0].split('_')
        if len(ref_info) == 1:
            genome_id = ref_info[0]
            species_real = ref_info[0]
        else:
            genome_id = '_'.join(ref_info[:3])
            species_real = '_'.join(ref_info[:2])
        if re.search('^K\d{1,}$', line_info[1]):
            gene_id = line_info[0]
            ko_num = line_info[1]
        elif re.search('^K\d{1,}$', line_info[2]):
            gene_id = line_info[1]
            ko_num = line_info[2]
        else:
            print (line)
            sys.exit()
        #gene_name = line_info[1]
        file_name = ''
        for file in glob.glob(f'{args.refdir}/*'):
            file_name = os.path.basename(file)
            if re.search(genome_id, file_name):
                genome_file = file
        if genome_file == '':
            continue
        if species_real not in species:
            print (species_real)
            continue
        seq_id, seq = getSeq(genome_file, gene_id)
        if seq_id == '' or seq == '':
            continue
        out_path = f'{args.outdir}/{ko_num}'
        os.makedirs(out_path) if not os.path.exists(out_path) else ''
        outfile = open(f'{out_path}/{os.path.basename(genome_file)}', 'a')
        #print (seq_id, genome_file, gene_id, seq)
        outfile.write(f'{seq_id} [{ko_num}] [{species[species_real]}]\n{seq}\n')
        outfile.close()
    analysis_out.close()
