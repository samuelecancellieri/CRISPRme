#crisprme download and test

echo "download all necessary folders and data"
#first download all the necessary data and test folders
wget https://github.com/samuelecancellieri/CRISPRme/raw/main/crisprme_test.tar.gz
tar -xzf crisprme_test.tar.gz

echo "unzip data"
#unzip annotations
cd crisprme_test/Annotations/
tar -xzf gencode_encode.hg38.tar.gz
cd ../

#unzip gencode
cd crisprme_test/Gencode/
tar -xzf gencode.protein_coding.tar.gz
cd ../

echo "start download VCF data and genome (this may take a long time due to connection speed)"
#download VCFs data
cd VCFs/
mkdir hg38_1000G
mkdir hg38_HGDP
#download 1000G
cd hg38_1000G/
echo "download 1000G VCFs"
for i in {1..22}
do
   wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr$i.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
done
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

#download HGDP
cd ../hg38_HGDP
echo "download HGDP VCFs"
for i in {1..22}
do
   wget ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chr$i.vcf.gz
done
wget ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chrX.vcf.gz
wget ftp://ngs.sanger.ac.uk:21/production/hgdp/hgdp_wgs.20190516/hgdp_wgs.20190516.full.chrY.vcf.gz
cd ../../

#download hg38
cd Genomes/
echo "download hg38"
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz
tar -xzf hg38.chromFa.tar.gz
mv chroms hg38

echo "start testing"
