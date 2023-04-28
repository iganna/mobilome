#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import pyreadr
from tqdm import tqdm
from os import listdir
from os.path import isfile, join
import numpy as np
import argparse
from argparse import ArgumentParser

parser = argparse.ArgumentParser()
parser.add_argument('--chr_thaliana', help='path to chr from annograma output')
parser.add_argument('--positions', help='path to chr corresponding positions in TAIR10')
parser.add_argument('--chr_lyrata', help='path to chr lyrata from annograma output')
parser.add_argument('--output_fixed_divergence',help='path to to the output of fixed divergence vcf')
parser.add_argument('--output_thaliana_converted_files',help='path to to the output thaliana snps vcf')
parser.add_argument('--output_thaliana_an_lyrata_converted_files',help='path to to the output thaliana and lyrata vcf')
parser.add_argument('--chr_number', help='Chromosome number in tair10')
args = parser.parse_args()

print('##### Reading files ######')
seqv_thaliana = pd.read_table(args.chr_thaliana, header=None, low_memory=False)
seqv_lyrata = pd.read_table(args.chr_lyrata, low_memory=False)
positions = pyreadr.read_r(args.positions)
print('##### DONE reading files ######')
print('##### Making alleles upper case ######')
# Make all letters uppercase ##IMPORTANT!!
seqv_thaliana = pd.concat([seqv_thaliana[col].astype(str).str.upper() for col in seqv_thaliana.columns], axis=1)
seqv_lyrata = pd.concat([seqv_lyrata[col].astype(str).str.upper()for col in seqv_lyrata.columns], axis=1)
print('##### Alleles converted to uppercase! ######')
# Extract positions as values (from readRDS function)
positions = list(positions.values())[0]
positions = positions.iloc[:, :-2]  # Remove unnesesary columns from positions
positions_first_col = positions['0']  # Exrtact positions of tair10
# Rename sequences columns as position columns
seqv_thaliana.columns = positions.columns

print('##### Merging positions ######')
# Merge positions and sequecnes columns
merged_df = pd.concat([positions_first_col.reset_index(drop=True), seqv_thaliana], axis=1)
merged_df.columns = ["positions", *merged_df.columns[1::]]
# Get dataframe without missing positions in tair10
without_missings_df = merged_df[merged_df['0'] != 'NAN']
#only_na_df = merged_df[~ merged_df['0'].notna()]
# only_na_df
# Make position integer
without_missings_df['positions'] = without_missings_df['positions'].astype('int')
seqv_lyrata['POS'] = seqv_lyrata['POS'].astype('int')
# Merge lyrata and thaliana
thaliana_lyrata_merged = without_missings_df.merge(seqv_lyrata, left_on='positions', right_on='POS')
# Drop POS column
thaliana_lyrata_merged = thaliana_lyrata_merged.drop(columns=['POS'])

#print('##### Calculating fixed divergence ######')

##FIXED DIVERGENCE###
#Iterate over columns of lyrata and thaliana and count differences
set_el_ar = 0
set_el_lyr = 0
el_thal = 0
difference_thal = []
for num, row in tqdm(enumerate(thaliana_lyrata_merged.itertuples(index=False)), total=thaliana_lyrata_merged.shape[0]):
    # for row in thaliana_lyrata_merged.itertuples(index=False):
    elements_arabidopsis = row[2:29]
    elements_lyrata = row[29:32]
    set_el_ar = set(elements_arabidopsis)
    set_el_lyr = set(elements_lyrata)
    if len(set_el_ar) != 1 or len(set_el_lyr) != 1:
        el_thal = '0'
    elif len(set_el_ar) == 1 and len(set_el_lyr) == 1:
        if set_el_lyr == {'NAN'} or set_el_ar == {'NAN'}:
            el_thal = '0'
        else:
            c = set_el_lyr.symmetric_difference(set_el_ar)
            if len(c) == 2:
                el_thal = '1'
    difference_thal.append(el_thal)

print('##### Calculating VCF thaliana+lyrata ######')
#Iterate  over rows and change reference allele to 0/0, alternative to 1/1 and missing to ./.
##lyrata and thaliana
##Taking into account multiallelic!!!!
chrom = args.chr_number  #Indicate chromosome number
resulting_path = args.output_thaliana_an_lyrata_converted_files
resulting_rows = []

for num, row in tqdm(enumerate(thaliana_lyrata_merged.itertuples()), total=thaliana_lyrata_merged.shape[0]):

    position = int(row[1])
    reference = row[2]
    vcf_id = num + 1
    alt = []
    out = []
    # A A A A A T NaN T -> A,T | expected: T
    #REF A, ALT T,G | expected 1/1,2/2
    for i in row[3:]:
        i=str(i)
        if i == "NAN":
            out.append("./.")
            continue
        elif i == reference:
            out.append("0/0")
        elif i!=reference:
            if i not in alt:
                alt.append(i)
                out.append(f'{alt.index(i)+1}/{alt.index(i)+1}')
            else:
                out.append(f'{alt.index(i)+1}/{alt.index(i)+1}')
    if len(alt)==0:
        alt.append('.')
    resulting_row = [chrom, position, vcf_id, reference, ",".join(list(alt)), ".", "PASS", ".", "GT", *out]
    resulting_rows.append(resulting_row)

print('##### Writing an output VCF thaliana+lyrata ######')
# Write resulting dataframe (positions from lyrata and thaliana)
resulting_df = pd.DataFrame(resulting_rows, columns=[*"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT".split(),
                                                     *thaliana_lyrata_merged.columns[2:]])
with open(f"{resulting_path}{chrom}_arabidopsis_processed_all_nucl_plus_lyr.vcf", "w") as out_file:
    out_file.write("##fileformat=VCFv4.0\n")
resulting_df.to_csv(f"{resulting_path}{chrom}_arabidopsis_processed_all_nucl_plus_lyr.vcf", sep="\t", index=False, mode="a")

# # Add fixed divergence in column
resulting_df['difference_thal'] = difference_thal
# # Extract only sites with fixed divergence (1)
resulting_df = resulting_df.loc[(resulting_df['difference_thal'] == '1')]
# # Drop unnecessary columns where we counted differences
resulting_df = resulting_df.drop(['difference_thal'], axis=1)

print('##### Writing an output Fixed divergence ######')

# Write only fixed diverged samples
resulting_path = args.output_fixed_divergence
with open(f"{resulting_path}{chrom}_arabidopsis_processed_all_nucl_plus_lyr_fixed_div.vcf", "w") as out_file:
    out_file.write("##fileformat=VCFv4.0\n")
resulting_df.to_csv(f"{resulting_path}{chrom}_arabidopsis_processed_all_nucl_plus_lyr_fixed_div.vcf", sep="\t", index=False, mode="a")

print('##### Calculating VCF only thaliana ######')
## Iterate  over rows and change reference allele to 0/0, alternative to 1/1 and missing to ./.
## Thaliana only data
chrom = args.chr_number  # Indicate chromosome number
resulting_path = args.output_thaliana_converted_files

resulting_rows = []

for num, row in tqdm(enumerate(without_missings_df.itertuples()), total=without_missings_df.shape[0]):
    position = int(row[1])
    reference = row[2]
    vcf_id = num + 1
    alt = []
    out = []
    # A A A A T NaN T -> A,T | expected: T
    #REF A, ALT T,G | expected 1/1,2/2
    for i in row[3:]:
        i=str(i)
        if i == "NAN":
            out.append("./.")
            continue
        elif i == reference:
            out.append("0/0")
        elif i!=reference:
            if i not in alt:
                alt.append(i)
                out.append(f'{alt.index(i)+1}/{alt.index(i)+1}')
            else:
                out.append(f'{alt.index(i)+1}/{alt.index(i)+1}')
    if len(alt)==0:
        alt.append('.')
    resulting_row = [chrom, position, vcf_id, reference, ",".join(list(alt)), ".", "PASS", ".", "GT", *out]
    resulting_rows.append(resulting_row)

print('##### Writing VCF only thaliana ######')
## Write resulting dataframe (only thaliana)
resulting_df = pd.DataFrame(resulting_rows, columns=[*"#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT".split(), *without_missings_df.columns[2:]])
with open(f"{resulting_path}{chrom}_arabidopsis_processed_all_nucl.vcf", "w") as out_file:
     out_file.write("##fileformat=VCFv4.0\n")
resulting_df.to_csv(f"{resulting_path}{chrom}_arabidopsis_processed_all_nucl.vcf", sep="\t", index=False, mode="a")
