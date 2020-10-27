#!/usr/bin/env python

import os 
import sys

vcf_dir = sys.argv[1]
files = os.listdir(vcf_dir)
vcf_name = sys.argv[2]
output_folder = sys.argv[3]

for f in files:
	print("Creating dictionary for file", f)
	number = f.split(".")[1].replace('chr', '')
	os.system("./creazione_dizionari_zipped.py "+vcf_dir+"/"+f+" "+str(number))
	os.system("mv "+'my_dict_chr' + number + '.json'+" "+output_folder+"/dictionaries_"+vcf_name+"/") 
