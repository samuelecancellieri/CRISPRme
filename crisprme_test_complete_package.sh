#!/bin/bash

echo "download the complete crisprme package and then unzip it"
axel -a -n 4 https://www.dropbox.com/s/pzfeb1k9v9ekyhr/complete_test_package_no_VCFs.tar.gz?dl=1
axel -a -n 4 https://www.dropbox.com/s/v1fxyhopjek6ib1/VCFs.tar.gz?dl=1
mkdir crisprme_complete_test
tar -xvzf complete_test_package_no_VCFs.tar.gz --directory crisprme_complete_test/
tar -xvzf VCFs.tar.gz --directory crisprme_complete_test/
cd crisprme_complete_test/

echo "generate input test data"
echo "CTAACAGTTGCTTTTATCACNNN" >sg1617.txt
echo "hg38_1000G/" >list_vcf.txt
echo "hg38_1000G.samplesID.txt" >list_samplesID.txt
echo "crisprme.py complete-search --genome Genomes/hg38/ --vcf list_vcf.txt/ --guide sg1617.txt --pam PAMs/20bp-NNN-SpCas9.txt --annotation Annotations/encode+gencode.hg38.bed --samplesID list_samplesID.txt --gene_annotation Annotations/gencode.protein_coding.bed --bMax 2 --mm 6 --bDNA 2 --bRNA 2 --merge 3 --output sg1617.6.2.2 --thread 8" >crisprme_complete_auto_test.sh

echo "to test the software, launch  crisprme_complete_auto_test.sh"
