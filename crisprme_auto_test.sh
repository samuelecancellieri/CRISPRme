#crisprme download and test

#first download all the necessary data and test folders
wget https://github.com/samuelecancellieri/CRISPRme/raw/main/crisprme_folders.zip
unzip crisprme_folders.zip

#unzip annotations
cd crisprme_test/Annotations/
unzip gencode_encode.hg38.zip
cd ../

#unzip gencode
cd crisprme_test/Gencode/
unzip gencode.protein_coding.zip
cd ../

#download VCFs data
cd VCFs/
mkdir hg38_1000G
mkdir hg38_HGDP
#download 1000G
cd hg38_1000G/
for i in {1..22}
do
   wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
done
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

#download HGDP
cd ../hg38_HGDP
for i in {1..22}
do
   wget ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr$i.vcf.gz
done
wget ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chrX.vcf.gz
wget ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chrY.vcf.gz
cd ../../

#download hg38
cd Genomes/
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar -xf hg38.chromFa.tar.gz
