#!/usr/bin/env python

from operator import truediv
from posixpath import expanduser
from intervaltree import IntervalTree
import sys
import time
import glob
import subprocess


def rev_comp(a):
    # perform reverse complement on a string
    if a == 'A' or a == 'a':
        return 'T'
    if a == 'T' or a == 't':
        return 'A'
    if a == 'C' or a == 'c':
        return 'G'
    return 'C'


def createBedforMultiAlternative(variantList, samples):
    # 5096 conta degli alleli di 1000G
    # bcftools view -H -r 1:222531770 -s HG00001 originalVCFs/ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
    # list containing snp position to search in the vcf
    toSearchList = list()
    # bed file with vcf data extracted
    vcfBedFile = outputDir + originFileName + '.temp.bed'
    # cycle over variant in list to extract the correct ones from vcf
    for elem in variantList:
        split = elem.strip().split('_')
        chrom = split[0]
        pos = split[1]
        vcfFile = ''
        for vcf in vcfList:
            if chrom+'.' in vcf:
                vcfFile = vcf
        vcfChr = chrom.replace('chr', '')
        toSearchList.append(vcfChr+':'+pos)
    # string with all the snps to search in the vcf file
    toSearchString = ','.join(toSearchList).strip()
    # subprocces run to call bcftools and extract the desidered snps
    subprocess.run(['bcftools', 'view', '-H', '-r',
                    toSearchString, '-s', samples, vcfFile, '-o', str(vcfBedFile)])
    haplotypeDict = dict()
    allele1 = list()
    allele2 = list()
    haplotype = 0
    # read line from vcf exctraction to obtain the alleles
    for line in open(vcfBedFile, 'r'):
        if 'INDEL' in line:
            continue
        split = line.strip().split('\t')
        tempHaploList = list()
        for elem in split:
            if '|' in elem:
                tempHaploList.append(elem)
        haplotypeDict[split[1]] = tempHaploList
        allele1 = [1]*len(tempHaploList)
        allele2 = [1]*len(tempHaploList)
    # count the correct haplotype frequence by and with sum operation
    # e.g. 0|1 and 1|1 will return 0|1 that will be counted as 1 for that sample haplotype
    for snp in haplotypeDict:
        for count, sample in enumerate(haplotypeDict[snp]):
            allele1_snp = sample.strip().split('|')[0]
            allele2_snp = sample.strip().split('|')[1]
            allele1[count] = int(allele1[count]) and int(allele1_snp)
            allele2[count] = int(allele2[count]) and int(allele2_snp)
    # final sum to generete haplotype frequency of the specific target
    for count, allele in enumerate(allele1):
        haplotype += int(allele1[count]) + int(allele2[count])
    # return to save the haplo freq divided by the number of alleles in 1000G (i.e. 2548samples*2alleles)
    return float(haplotype)/5096


print('READING INPUT FILES')

# READ INPUT FILES
crispritzResultFile = sys.argv[1]  # file with CRISPRitz results
empiricalResults = sys.argv[2]  # file with empirical results
annotationFile = sys.argv[3]  # annotation file used to find nearest gene
guideFile = sys.argv[4]  # real guide used in the search
outputDir = sys.argv[5]  # output file name
check = sys.argv[6].upper()  # check if user input a GENCODE annotation file
genomeRelease = str(sys.argv[7]).strip()  # genome used in the search phase
# directory of vcf to perform haplotype count with more than one SVs in single target
vcfFileDirectory = sys.argv[8]

# OPEN INPUT FILES AND PREPARE OUTPUT FILE
# crispritz results file open
inCrispritzResults = open(crispritzResultFile, 'r')
# empirical-seq data open
inEmpiricalResults = open(empiricalResults, 'r').readlines()
# annotation data open
inAnnotationFile = open(annotationFile, 'r')
# guide file
# realguide = open(guideFile, 'r').readlines()
# open outputDir to write results
originFileName = crispritzResultFile.split('/')
originFileName = originFileName[len(originFileName)-1]
outFile = open(outputDir + originFileName +
               '.integrated_results.tsv', 'w')
# list file in vcf directory
checkVCF = False
if 'vuota' not in vcfFileDirectory:
    checkVCF = True
    vcfList = glob.glob(vcfFileDirectory+'/*.gz')
    redirectFile = open(outputDir + originFileName + '.redirectFile.out', 'w')
    for vcfFile in vcfList:
        # create tabix index if not already done
        subprocess.run(['tabix', str(vcfFile)], stderr=redirectFile)
    redirectFile.close()


empiricalTree = IntervalTree()
empiricalList = []
genomeDict = {}
empiricalDict = {}
valueDict = {}
# saveDict = {"real_guide": 'n', "genome": 'n', "chr": 'n', "prim_pos": 'n', "strand": 'n', "highest_CFD_guide_alignment": 'n', "highest_CFD_alignment(ref)": 'n',
#             "highest_CFD_alignment(alt)": 'n', "ref_seq_length": 'n', "ref_pos_alt(aligned_strand)": 'n', "pam": 'n', "annotation": 'n', "highest_CFD_score": 'n',
#             "highest_CFD_score(ref)": 'n', "highest_CFD_score(alt)": 'n', "risk_score": 'n', "absolute_risk_score": 'n', "highest_CFD_mismatch": 'n',
#             "highest_CFD_bulge": 'n', "highest_CFD_mismatch+bulge": 'n', "fewest_mm+bulge_guide_alignment": 'n', "fewest_mm+bulge_alignment(ref)": 'n',
#             "fewest_mm+bulge_alignment(alt)": 'n', "fewest_mm+bulge_CFD_score(ref)": 'n', "fewest_mm+bulge_CFD_score(alt)": 'n', "fewest_mismatch": 'n',
#             "fewest_bulge": 'n', "fewest_mismatch+bulge": 'n', "alt_haplotypes": 'n', "prim_origin": 'n', "prim_AF": 'n', "prim_samples": 'n',
#             "prim_SNP_ID(positive_strand)": 'n', "gene_name": 'n', "gene_ID": 'n', "gene_annotation": 'n', "gene_distance(kb)": 'n', "lowest_empirical": 'n',
#             "Nature2019": 'n', "Nature2019_mm+bul": 'n', "CHANGEseq": 'n', "CHANGEseq_mm+bul": 'n', "CIRCLEseq": 'n', "CIRCLEseq_mm+bul": 'n', "ONEseq": 'n',
#             "ONEseq_mm+bul": 'n', "GUIDEseq_293": 'n', "GUIDEseq_293_mm+bul": 'n', "GUIDEseq_CD34": 'n', "GUIDEseq_CD34_mm+bul": 'n', "GUIDEseq": 'n',
#             "GUIDEseq_mm+bul": 'n'}

# saveDict = {"real_guide": 'n', "genome": 'n', "chr": 'n', "prim_pos": 'n', "strand": 'n', "highest_CFD_guide_alignment": 'n', "highest_CFD_alignment(ref)": 'n',
#             "highest_CFD_alignment(alt)": 'n', "ref_seq_length": 'n', "ref_pos_alt(aligned_strand)": 'n', "pam": 'n', "annotation": 'n', "highest_CFD_score": 'n',
#             "highest_CFD_score(ref)": 'n', "highest_CFD_score(alt)": 'n', "risk_score": 'n', "absolute_risk_score": 'n', "highest_CFD_mismatch": 'n',
#             "highest_CFD_bulge": 'n', "highest_CFD_mismatch+bulge": 'n', "fewest_mm+bulge_guide_alignment": 'n', "fewest_mm+bulge_alignment(ref)": 'n',
#             "fewest_mm+bulge_alignment(alt)": 'n', "fewest_mm+bulge_CFD_score(ref)": 'n', "fewest_mm+bulge_CFD_score(alt)": 'n', "fewest_mismatch": 'n',
#             "fewest_bulge": 'n', "fewest_mismatch+bulge": 'n', "alt_haplotypes": 'n', "prim_origin": 'n', "prim_AF": 'n', "prim_samples": 'n',
#             "prim_SNP_ID(positive_strand)": 'n', "gene_name": 'n', "gene_ID": 'n', "gene_annotation": 'n', "gene_distance(kb)": 'n', "lowest_empirical": 'n'}

saveDict = {
    'Highest_CFD_Strand': 'NA',
    'Chromosome': 'NA',
    'Highest_CFD_start_coordinate': 'NA',
    'Highest_CFD_aligned_spacer+PAM': 'NA',
    'Highest_CFD_aligned_protospacer+PAM_REF': 'NA',
    'Highest_CFD_aligned_protospacer+PAM_ALT': 'NA',
    'Highest_CFD_PAM': 'NA',
    'Highest_CFD_mismatches': 'NA',
    'Highest_CFD_bulges': 'NA',
    'Highest_CFD_mismatches+bulges': 'NA',
    'Highest_CFD_bulge_type': 'NA',
    'Highest_CFD_REF/ALT_origin': 'NA',
    'Highest_CFD_PAM_gen': 'NA',
    'Highest_CFD_score': 'NA',
    'Highest_CFD_score_REF': 'NA',
    'Highest_CFD_score_ALT': 'NA',
    'Highest_CFD_risk_score': 'NA',
    'Highest_CFD_variant_info_spacer+PAM': 'NA',
    'Highest_CFD_variant_info_genome': 'NA',
    'Highest_CFD_variant_MAF': 'NA',
    'Highest_CFD_variant_rsID': 'NA',
    'Highest_CFD_variant_samples': 'NA',
    'Not_found_in_REF': 'NA',
    'Other_motifs': 'NA',
    'Fewest_mm+b_Strand': 'NA',
    'Fewest_mm+b_start_coordinate': 'NA',
    'Fewest_mm+b_start_coordinate': 'NA',
    'Fewest_mm+b_aligned_spacer+PAM': 'NA',
    'Fewest_mm+b_aligned_protospacer+PAM_REF': 'NA',
    'Fewest_mm+b_aligned_protospacer+PAM_ALT': 'NA',
    'Fewest_mm+b_PAM': 'NA',
    'Fewest_mm+b_mismatches': 'NA',
    'Fewest_mm+b_bulges': 'NA',
    'Fewest_mm+b_mismatches+bulges': 'NA',
    'Fewest_mm+b_bulge_type': 'NA',
    'Fewest_mm+b_REF/ALT_origin': 'NA',
    'Fewest_mm+b_PAM_gen': 'NA',
    'Fewest_mm+b_CFD_score': 'NA',
    'Fewest_mm+b_CFD_score_REF': 'NA',
    'Fewest_mm+b_CFD_score_ALT': 'NA',
    'Fewest_mm+b_CFD_risk_score': 'NA',
    'Fewest_mm+b_variant_info_spacer+PAM': 'NA',
    'Fewest_mm+b_variant_info_genome': 'NA',
    'Fewest_mm+b_variant_MAF': 'NA',
    'Fewest_mm+b_variant_rsID': 'NA',
    'Fewest_mm+b_variant_samples': 'NA',
    'Annotation_GENCODE': 'NA',
    'Annotation_closest_gene_name': 'NA',
    'Annotation_closest_gene_ID': 'NA',
    'Annotation_closest_gene_distance_(kb)': 'NA',
    'Annotation_ENCODE': 'NA',
    'Annotation_personal': 'NA'
}


start_time = time.time()

print('CREATING INTERVAL TREES')

for count, line in enumerate(inEmpiricalResults):
    empList = line.strip().split('\t')
    empList = [elem.strip() for elem in empList]
    empList.append(count)
    # adding empirical data to the tree
    empiricalTree[int(empList[1]):int(empList[2])] = empList
    # generate temp empirical dict to save empirical information per row
    empiricalDict[str(empList[4])] = 50
    # value of empirical row, keeping info about mm+bul for each empirical origin (seq data, in-silico, vitro, ecc)
    valueDict[str(empList[5])] = 'NA'
    # update save dict with user-defined names from empirical data
    saveDict[str(empList[5])] = 'NA'
    # newkey = str(empList[5])+'_mm+bul'
    # saveDict[newkey] = 'NA'

# writing header in file
save = '#'
save += '\t'.join(list(saveDict.keys()))
# for key in saveDict:
#     save += str(key)+'\t'
save += '\n'
outFile.write(save)

print('INTEGRATING RESULTS')

if '#' in inCrispritzResults.readline():
    print('SKIP HEADER')
else:
    inCrispritzResults.seek(0)

for nline, line in enumerate(inCrispritzResults):
    x = line.strip().split('\t')
    try:
        annotationLine = inAnnotationFile.readline().strip().split('\t')
    except:
        annotationFile = 'NA'
    # print(annotationLine)
    lowestEmpirical = 100

    for key in saveDict:
        saveDict[key] = 'NA'

    for key in empiricalDict:
        empiricalDict[key] = 50
        valueDict[key] = 'NA'

    if 'NA' not in annotationLine and check == 'TRUE':
        for elem in annotationLine:
            if 'gene_id' in elem:
                temp = elem.strip().split(';')
                for name in temp:
                    if 'gene_id' in name:
                        saveDict['Annotation_closest_gene_ID'] = name.strip().split('=')[
                            1]
                    if 'gene_name' in name:
                        saveDict['Annotation_closest_gene_name'] = name.strip().split('=')[
                            1]
                    # if 'gene_type' in name:
                    #     tipo = str(annotationLine[11])
                    #     if 'codon' in tipo or 'exon' in tipo:
                    #         tipo = 'CDS'
                    #     elif 'gene' in tipo or 'transcript' in tipo:
                    #         tipo = 'intron'
                    #     elif 'five_prime_UTR' in tipo:
                    #         tipo = "5'UTR"
                    #     elif 'three_prime_UTR' in tipo:
                    #         tipo = "3'UTR"
                    # saveDict['gene_annotation'] = str(annotationLine[11])
                        saveDict['Annotation_GENCODE'] = str(
                            annotationLine[11])
        saveDict['Annotation_closest_gene_distance_(kb)'] = str(
            float(annotationLine[len(annotationLine)-1])/1000)
        if float(annotationLine[len(annotationLine)-1]) != 0:
            saveDict['Annotation_GENCODE'] = 'intergenic'

    # origin = ''
    # if 'n' in str(x[13]):
    #     origin = 'ref'
    # else:
    #     origin = 'alt'

    variantList = ['NA']
    if str(x[13]) != 'NA' and str(x[13]) != 'n':
        variantList = str(x[18]).strip().split(',')
        if len(variantList) > 1 and checkVCF:
            samples = x[13]
            x[17] = createBedforMultiAlternative(variantList, samples)
            x[17] = str(x[17])[:7]
        var_pos = []
        # generate variant position corrected to be in the positive strand
        if '+' in str(x[7]):
            refseq = str(x[3])
            altseq = str(x[2])
            for pos, nucleotide in enumerate(refseq):
                if altseq[pos].lower() != nucleotide.lower():
                    var_pos.append(pos+1)
        else:
            refseq = str(x[3])
            altseq = str(x[2])
            for pos in range(len(refseq)-1, -1, -1):
                if refseq[pos].lower() != altseq[pos].lower():
                    var_pos.append(pos+1)
        # generate variant nucleotide reverse complemented always in positive strand
        for count, elem in enumerate(variantList):
            split = str(elem).strip().split('_')
            split_one_len = len(split[2])
            split_second_len = len(split[3])
            correction = 0
            if split_one_len != split_second_len:
                correction = split_second_len
            if '+' in str(x[7]):
                variantList[count] = str(
                    split[2])+str(int(var_pos[count])+correction)+str(split[3])
            else:
                firstcomp = ''
                secondcomp = ''
                for piece in str(split[2]):
                    firstcomp += rev_comp(piece)
                for piece in str(split[3]):
                    secondcomp += rev_comp(piece)
                variantList[count] = ''.join(
                    reversed(firstcomp))+str(int(var_pos[count])+correction)+''.join(reversed(secondcomp))
    variantList_highest_cfd = variantList

    variantList = ['NA']
    if str(x[37]) != 'NA' and str(x[37]) != 'n':
        variantList = str(x[42]).strip().split(',')
        if len(variantList) > 1 and checkVCF:
            samples = x[37]
            x[41] = createBedforMultiAlternative(variantList, samples)
            x[41] = str(x[41])[:7]
        var_pos = []
        # generate variant position corrected to be in the positive strand
        if '+' in str(x[31]):
            refseq = str(x[27])
            altseq = str(x[26])
            for pos, nucleotide in enumerate(refseq):
                if altseq[pos].lower() != nucleotide.lower():
                    var_pos.append(pos+1)
        else:
            refseq = str(x[27])
            altseq = str(x[26])
            for pos in range(len(refseq)-1, -1, -1):
                if refseq[pos].lower() != altseq[pos].lower():
                    var_pos.append(pos+1)
        # generate variant nucleotide reverse complemented always in positive strand
        for count, elem in enumerate(variantList):
            split = str(elem).strip().split('_')
            split_one_len = len(split[2])
            split_second_len = len(split[3])
            correction = 0
            if split_one_len != split_second_len:
                correction = split_second_len
            if '+' in str(x[7]):
                variantList[count] = str(
                    split[2])+str(int(var_pos[count])+correction)+str(split[3])
            else:
                firstcomp = ''
                secondcomp = ''
                for piece in str(split[2]):
                    firstcomp += rev_comp(piece)
                for piece in str(split[3]):
                    secondcomp += rev_comp(piece)
                variantList[count] = ''.join(
                    reversed(firstcomp))+str(int(var_pos[count])+correction)+''.join(reversed(secondcomp))
    variantList_fewest_mm_b = variantList

    saveDict['Highest_CFD_Strand'] = str(x[7])
    saveDict['Chromosome'] = str(x[4])
    saveDict['Highest_CFD_start_coordinate'] = str(x[5])
    saveDict['Highest_CFD_aligned_spacer+PAM'] = str(x[1])
    saveDict['Highest_CFD_aligned_protospacer+PAM_REF'] = str(x[3])
    saveDict['Highest_CFD_aligned_protospacer+PAM_ALT'] = str(x[2])
    saveDict['Highest_CFD_PAM'] = 'NA'
    saveDict['Highest_CFD_mismatches'] = str(x[8])
    saveDict['Highest_CFD_bulges'] = str(x[9])
    saveDict['Highest_CFD_mismatches+bulges'] = str(x[10])
    saveDict['Highest_CFD_bulge_type'] = str(x[0])
    saveDict['Highest_CFD_REF/ALT_origin'] = 'ref' if str(
        x[13]) == 'NA' else 'alt'
    saveDict['Highest_CFD_PAM_gen'] = str(x[11])
    saveDict['Highest_CFD_score'] = str(x[20]) if float(
        x[20]) > float(x[21]) else str(x[21])
    saveDict['Highest_CFD_score_REF'] = str(x[21])
    saveDict['Highest_CFD_score_ALT'] = str(x[20])
    saveDict['Highest_CFD_risk_score'] = str(x[22])
    saveDict['Highest_CFD_variant_info_spacer+PAM'] = ','.join(
        variantList_highest_cfd)
    saveDict['Highest_CFD_variant_info_genome'] = str(x[18])
    saveDict['Highest_CFD_variant_MAF'] = str(x[17])
    saveDict['Highest_CFD_variant_rsID'] = 'NA' if str(
        x[16]) == '.' else str(x[16])
    saveDict['Highest_CFD_variant_samples'] = str(x[13])
    saveDict['Not_found_in_REF'] = 'y' if str(x[12]) == 'y' else 'NA'
    saveDict['Other_motifs'] = str(x[19])
    saveDict['Fewest_mm+b_Strand'] = str(x[31])
    saveDict['Fewest_mm+b_start_coordinate'] = str(x[29])
    saveDict['Fewest_mm+b_aligned_spacer+PAM'] = str(x[25])
    saveDict['Fewest_mm+b_aligned_protospacer+PAM_REF'] = str(x[27])
    saveDict['Fewest_mm+b_aligned_protospacer+PAM_ALT'] = str(x[26])
    saveDict['Fewest_mm+b_PAM'] = 'NA'
    saveDict['Fewest_mm+b_mismatches'] = str(x[32])
    saveDict['Fewest_mm+b_bulges'] = str(x[33])
    saveDict['Fewest_mm+b_mismatches+bulges'] = str(x[34])
    saveDict['Fewest_mm+b_bulge_type'] = str(x[24])
    saveDict['Fewest_mm+b_REF/ALT_origin'] = 'ref' if str(
        x[37]) == 'NA' else 'alt'
    saveDict['Fewest_mm+b_PAM_gen'] = str(x[35])
    saveDict['Fewest_mm+b_CFD_score'] = str(x[44]
                                            ) if float(x[44]) > float(x[45]) else str(x[45])
    saveDict['Fewest_mm+b_CFD_score_REF'] = str(x[45])
    saveDict['Fewest_mm+b_CFD_score_ALT'] = str(x[44])
    saveDict['Fewest_mm+b_CFD_risk_score'] = str(x[46])
    saveDict['Fewest_mm+b_variant_info_spacer+PAM'] = ','.join(
        variantList_fewest_mm_b)
    saveDict['Fewest_mm+b_variant_info_genome'] = str(x[42])
    saveDict['Fewest_mm+b_variant_MAF'] = str(x[41])
    saveDict['Fewest_mm+b_variant_rsID'] = 'NA' if str(
        x[40]) == '.' else str(x[40])
    saveDict['Fewest_mm+b_variant_samples'] = str(x[37])

    if saveDict['Highest_CFD_REF/ALT_origin'] == 'ref':
        saveDict['Highest_CFD_aligned_protospacer+PAM_REF'] = saveDict['Highest_CFD_aligned_protospacer+PAM_ALT']
        saveDict['Highest_CFD_aligned_protospacer+PAM_ALT'] = 'NA'
    if saveDict['Fewest_mm+b_REF/ALT_origin'] == 'ref':
        saveDict['Fewest_mm+b_aligned_protospacer+PAM_REF'] = saveDict['Fewest_mm+b_aligned_protospacer+PAM_ALT']
        saveDict['Fewest_mm+b_aligned_protospacer+PAM_ALT'] = 'NA'

    change_alt_ref_highest_cfd = False
    if saveDict['Highest_CFD_REF/ALT_origin'] == 'alt' and saveDict['Highest_CFD_score_REF'] == saveDict['Highest_CFD_score_ALT']:
        if 'DNA' in saveDict['Highest_CFD_bulge_type']:
            mm = 0
            bulge = 0
            for nt in saveDict['Highest_CFD_aligned_protospacer+PAM_REF']:
                if nt.islower():
                    mm += 1
            for nt in saveDict['Highest_CFD_aligned_spacer+PAM']:
                if nt == '-':
                    bulge += 1
        else:
            mm = 0
            bulge = 0
            for nt in saveDict['Highest_CFD_aligned_protospacer+PAM_REF']:
                if nt.islower():
                    mm += 1
                if nt == '-':
                    bulge += 1
        if mm <= int(saveDict['Highest_CFD_mismatches']) and bulge <= int(saveDict['Highest_CFD_bulges']):
            change_alt_ref_highest_cfd = True

    change_alt_ref_fewest_mm_b = False
    if saveDict['Fewest_mm+b_REF/ALT_origin'] == 'alt' and saveDict['Fewest_mm+b_CFD_score_REF'] == saveDict['Fewest_mm+b_CFD_score_ALT']:
        if 'DNA' in saveDict['Fewest_mm+b_bulge_type']:
            mm = 0
            bulge = 0
            for nt in saveDict['Fewest_mm+b_aligned_protospacer+PAM_REF']:
                if nt.islower():
                    mm += 1
            for nt in saveDict['Fewest_mm+b_aligned_spacer+PAM']:
                if nt == '-':
                    bulge += 1
        else:
            mm = 0
            bulge = 0
            for nt in saveDict['Fewest_mm+b_aligned_protospacer+PAM_REF']:
                if nt.islower():
                    mm += 1
                if nt == '-':
                    bulge += 1
        if mm <= int(saveDict['Fewest_mm+b_mismatches']) and bulge <= int(saveDict['Fewest_mm+b_bulges']):
            change_alt_ref_fewest_mm_b = True

    if change_alt_ref_highest_cfd:
        saveDict['Highest_CFD_aligned_protospacer+PAM_REF'] = saveDict['Highest_CFD_aligned_protospacer+PAM_ALT']
        saveDict['Highest_CFD_aligned_protospacer+PAM_ALT'] = 'NA'
        saveDict['Highest_CFD_REF/ALT_origin'] = 'ref'
        saveDict['Highest_CFD_variant_info_spacer+PAM'] = 'NA'
        saveDict['Highest_CFD_variant_info_genome'] = 'NA'
        saveDict['Highest_CFD_variant_MAF'] = 'NA'
        saveDict['Highest_CFD_variant_rsID'] = 'NA'
        saveDict['Highest_CFD_variant_samples'] = 'NA'

    if change_alt_ref_fewest_mm_b:
        saveDict['Fewest_mm+b_aligned_protospacer+PAM_REF'] = saveDict['Fewest_mm+b_aligned_protospacer+PAM_ALT']
        saveDict['Fewest_mm+b_aligned_protospacer+PAM_ALT'] = 'NA'
        saveDict['Fewest_mm+b_REF/ALT_origin'] = 'ref'
        saveDict['Fewest_mm+b_variant_info_spacer+PAM'] = 'NA'
        saveDict['Fewest_mm+b_variant_info_genome'] = 'NA'
        saveDict['Fewest_mm+b_variant_MAF'] = 'NA'
        saveDict['Fewest_mm+b_variant_rsID'] = 'NA'
        saveDict['Fewest_mm+b_variant_samples'] = 'NA'

    count_N_in_guide = 0  # check how long is the pam counting Ns in the guide
    pam_at_start = False  # check if pam is at start of the sequence
    # count number of Ns in the guide
    for count, elem in enumerate(str(x[15])):
        if elem == 'N':
            count_N_in_guide += 1
            if count == 0:  # if N is at start of the guide, pam_at_start = true
                pam_at_start = True

    if saveDict['Highest_CFD_REF/ALT_origin'] == 'ref':
        if pam_at_start:  # save pam sequence extracting directly from the ref sequence
            saveDict['Highest_CFD_PAM'] = saveDict['Highest_CFD_aligned_protospacer+PAM_REF'][:count_N_in_guide]
        else:
            saveDict['Highest_CFD_PAM'] = saveDict['Highest_CFD_aligned_protospacer+PAM_REF'][-count_N_in_guide:]
    else:
        if pam_at_start:  # save pam sequence extracting directly from the var sequence
            saveDict['Highest_CFD_PAM'] = saveDict['Highest_CFD_aligned_protospacer+PAM_ALT'][:count_N_in_guide]
        else:
            saveDict['Highest_CFD_PAM'] = saveDict['Highest_CFD_aligned_protospacer+PAM_ALT'][-count_N_in_guide:]

    if saveDict['Fewest_mm+b_REF/ALT_origin'] == 'ref':
        if pam_at_start:  # save pam sequence extracting directly from the ref sequence
            saveDict['Fewest_mm+b_PAM'] = saveDict['Fewest_mm+b_aligned_protospacer+PAM_REF'][:count_N_in_guide]
        else:
            saveDict['Fewest_mm+b_PAM'] = saveDict['Fewest_mm+b_aligned_protospacer+PAM_REF'][-count_N_in_guide:]
    else:
        if pam_at_start:  # save pam sequence extracting directly from the var sequence
            saveDict['Fewest_mm+b_PAM'] = saveDict['Fewest_mm+b_aligned_protospacer+PAM_ALT'][:count_N_in_guide]
        else:
            saveDict['Fewest_mm+b_PAM'] = saveDict['Fewest_mm+b_aligned_protospacer+PAM_ALT'][-count_N_in_guide:]

    annotationList = x[14].split(',')
    personal_annotations = list()
    encode_annotations = list()
    for elem in annotationList:
        if '_personal' not in elem:
            encode_annotations.append(elem)
        else:
            personal_annotations.append(elem)

    if len(personal_annotations) > 0:
        saveDict['Annotation_personal'] = ','.join(personal_annotations)
    if len(encode_annotations) > 0:
        saveDict['Annotation_ENCODE'] = ','.join(encode_annotations)

    foundEmpirical = sorted(empiricalTree[int(x[6])-4:int(x[6])+4])

    for found in range(0, len(foundEmpirical)):
        empirical = foundEmpirical[found].data
        if str(saveDict['chr']) == str(empirical[0]):
            empiricalList.append(empirical[7])
            valueDict[str(empirical[4])] = empirical[5]
            empiricalDict[str(empirical[4])] = int(empirical[3])

    for key in empiricalDict:
        if int(empiricalDict[key]) < 50:
            saveDict[key] = str(valueDict[key])
            newkey = str(key)+'_mm+bul'
            saveDict[newkey] = empiricalDict[key]
            # if int(empiricalDict[key]) < lowestEmpirical:
            # saveDict['lowest_empirical'] = str(empiricalDict[key])

    save = ''
    for key in saveDict:
        save += str(saveDict[key])+'\t'
    save += '\n'

    outFile.write(save)

print('CHECKING MISSING RESULTS')

notFoundFile = open(outputDir + originFileName +
                    '.empirical_not_found.tsv', 'w')

for count, line in enumerate(inEmpiricalResults):
    if count not in empiricalList:
        notFoundFile.write(line)

print("INTEGRATION COMPLETED IN: %s seconds" % (time.time() - start_time))
