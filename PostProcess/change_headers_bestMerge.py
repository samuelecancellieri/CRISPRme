#!/usr/bin/env python
"""
Created on Sat May 29 18:02:45 2021

@author: franc
"""

'''
Script used to convert from old bestMerge format to new alt_results format
'''
import pandas as pd
import numpy as np
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

in_path = sys.argv[1]
out_path = sys.argv[2]

chunksize_ = 5000000000000

# if not isFileLocked(in_path):
chunks = pd.read_csv(in_path, sep='\t', chunksize=chunksize_, na_filter=False)

new_order = ['Real_Guide', 'Direction', 'Chromosome', 'Position', 'Cluster_Position', 'crRNA', 'Reference', 'DNA', 'Mismatches', 'Bulge_Size', 'Total',
             '#Bulge_type', 'PAM_gen', 'CFD', 'CFD_ref', 'Highest_CFD_Risk_Score', 'Highest_CFD_Absolute_Risk_Score', 'Var_uniq', 'SNP', 'AF', 'rsID',
             'Samples', '#Seq_in_cluster', 'MMBLG_Real_Guide', 'MMBLG_Chromosome', 'MMBLG_Position', 'MMBLG_Cluster_Position', 'MMBLG_Direction', 'MMBLG_crRNA',
             'MMBLG_Reference', 'MMBLG_DNA', 'MMBLG_Mismatches', 'MMBLG_Bulge_Size', 'MMBLG_Total', 'MMBLG_#Bulge_type', 'MMBLG_PAM_gen', 'MMBLG_CFD', 'MMBLG_CFD_ref',
             'MMBLG_CFD_Risk_Score', 'MMBLG_CFD_Absolute_Risk_Score', 'MMBLG_Var_uniq', 'MMBLG_SNP', 'MMBLG_AF', 'MMBLG_rsID', 'MMBLG_Samples', 'MMBLG_#Seq_in_cluster',
             'MMBLG_Annotation_Type', 'Annotation_Type']

to_remove = ['Real_Guide', 'Cluster_Position', 'Highest_CFD_Absolute_Risk_Score', 'MMBLG_Real_Guide', 'MMBLG_Chromosome', 'MMBLG_Cluster_Position',
             'MMBLG_CFD_Absolute_Risk_Score', 'MMBLG_Var_uniq', 'MMBLG_#Seq_in_cluster', 'MMBLG_Annotation_Type'
             ]

new_names = ['Highest_CFD_Strand', 'Chromosome', 'Start_coordinate', 'Highest_CFD_aligned_spacer+PAM', 'Highest_CFD_aligned_protospacer+PAM_REF', 'Highest_CFD_aligned_protospacer+PAM_ALT',
             'Highest_CFD_mismatches', 'Highest_CFD_bulges', 'Highest_CFD_mismatches+bulges', 'Highest_CFD_bulge_type', 'Highest_CFD_PAM_gen', 'Highest_CFD_score', 'Highest_CFD_score_REF',
             'Highest_CFD_risk_score', 'Not_found_in_REF', 'Variant_info_genome', 'Highest_CFD_variant_MAF', 'Highest_CFD_variant_rsID', 'Highest_CFD_variant_Samples', 'Other_motifs', 'Fewest_mm+b_Strand', 'Fewest_mm+b_start_coordinate',
             'Fewest_mm+b_aligned_spacer+PAM', 'Fewest_mm+b_aligned_protospacer+PAM_REF', 'Fewest_mm+b_aligned_protospacer+PAM_ALT', 'Fewest_mm+b_mismatches', 'Fewest_mm+b_bulges',
             'Fewest_mm+b_mismatches+bulges', 'Fewest_mm+b_bulge_type', 'Fewest_mm+b_PAM_gen', 'Fewest_mm+b_CFD_score', 'Fewest_mm+b_CFD_score_REF', 'Fewest_mm+b_CFD_risk_score',
             'Fewest_mm+b_variant_info_genome', 'Fewest_mm+b_variant_MAF', 'Fewest_mm+b_variant_rsID', 'Fewest_mm+b_variant_samples', 'Annotation_ENCODE']


header = True
for chunk in chunks:

    chunk = chunk[new_order]
    chunk = chunk.drop(to_remove, axis=1)
    chunk.columns = new_names
    #chunk[chunk.columns.difference(['Not_found_in_REF'])] = chunk[chunk.columns.difference(['Not_found_in_REF'])].replace('n', 'NA')
    chunk = chunk.replace('n', 'NA')
    #chunk = chunk.replace(regex=['\*.,\*', '\*,.\*'], value='NA')
    chunk['rsID'] = chunk['rsID'].str.replace('.', 'NA')
    mask = chunk['Highest_CFD_aligned_protospacer+PAM_REF'] == 'NA'
    chunk['Highest_CFD_aligned_protospacer+PAM_REF_corrected'] = np.where(
        mask, chunk['Highest_CFD_aligned_protospacer+PAM_ALT'], chunk['Highest_CFD_aligned_protospacer+PAM_REF'])
    chunk['Highest_CFD_aligned_protospacer+PAM_ALT'] = np.where(
        mask, chunk['Highest_CFD_aligned_protospacer+PAM_REF'], chunk['Highest_CFD_aligned_protospacer+PAM_ALT'])

    chunk['Highest_CFD_aligned_protospacer+PAM_REF'] = chunk['Highest_CFD_aligned_protospacer+PAM_REF_corrected']
    chunk.drop('Highest_CFD_aligned_protospacer+PAM_REF_corrected',
               axis=1, inplace=True)

    mask = chunk['Fewest_mm+b_aligned_protospacer+PAM_REF'] == 'NA'
    chunk['Fewest_mm+b_aligned_protospacer+PAM_REF_corrected'] = np.where(
        mask, chunk['Fewest_mm+b_aligned_protospacer+PAM_ALT'], chunk['Fewest_mm+b_aligned_protospacer+PAM_REF'])
    chunk['Fewest_mm+b_aligned_protospacer+PAM_ALT'] = np.where(
        mask, chunk['Fewest_mm+b_aligned_protospacer+PAM_REF'], chunk['Fewest_mm+b_aligned_protospacer+PAM_ALT'])

    chunk['Fewest_mm+b_aligned_protospacer+PAM_REF'] = chunk['Fewest_mm+b_aligned_protospacer+PAM_REF_corrected']
    chunk.drop('Fewest_mm+b_aligned_protospacer+PAM_REF_corrected',
               axis=1, inplace=True)

    chunk.to_csv(out_path, header=header, mode='w',
                 sep='\t', index=False, na_rep='NA')

    header = False
