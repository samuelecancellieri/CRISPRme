#!/usr/bin/env python

import os
import sys
from multiprocessing import Pool

ref_folder = sys.argv[1]
vcf_dir = sys.argv[2]
vcf_name = sys.argv[3]
guide_file = sys.argv[4]
pam_file = sys.argv[5]
bMax = sys.argv[6]
mm = sys.argv[7]
bDNA = sys.argv[8]
bRNA = sys.argv[9]
output_folder = sys.argv[10]

global use_thread 
use_thread = 0

def start_indels(f):
    global use_thread
    splitted = f.split('.')
    chrom = splitted[1]
    if use_thread == t-1:
        use_thread = 0
    os.system(f"./indels_process.sh \"{vcf_dir}/{f}\" \"{chrom}\" \"{output_folder}\" \"{vcf_name}\" \"{ref_folder}\" \"{pam_file}\" \"{guide_file}\" {bMax} {mm} {bDNA} {bRNA} {use_thread}")
    use_thread =+ 1

chrs = []
for f in os.listdir(vcf_dir):
    if 'vcf.gz' in f:
        chrs.append(f)

cpus = len(os.sched_getaffinity(0))
if cpus - 3 < 10:
    if cpus - 3 < 0:
        t = 1
    else:
        t = cpus - 3
else:
    t = 10

with Pool(processes = t) as pool:
    pool.map(start_indels, chrs)