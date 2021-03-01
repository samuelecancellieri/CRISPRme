#!/usr/bin
import sys
import json

inFasta = open(sys.argv[1], 'r')  # lettura fasta del chr
chr = inFasta.readline().strip().replace('>', '')  # lettura fasta del chr
genomeStr = inFasta.readlines()  # lettura fasta del chr
genomeStr = ''.join(genomeStr).upper()
genomeStr = genomeStr.replace('\n', '')
inDict = open(sys.argv[2], 'r')
inTarget = open(sys.argv[3], 'r')

mydict = json.load(inDict)


iupac_code = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "r": ("A", "G"),
    "y": ("C", "T"),
    "s": ("G", "C"),
    "w": ("A", "T"),
    "k": ("G", "T"),
    "m": ("A", "C"),
    "b": ("C", "G", "T"),
    "d": ("A", "G", "T"),
    "h": ("A", "C", "T"),
    "v": ("A", "C", "G"),
    'N': ('A', 'T', 'C', 'G')
}

# For scoring of CFD And Doench
tab = str.maketrans("ACTGRYSWMKHDBVactgryswmkhdbv",
                    "TGACYRSWKMDHVBtgacyrswkmdhvb")


def reverse_complement_table(seq):
    return seq.translate(tab)[::-1]


def retrieveFromDict(chr_pos):
    entry = mydict[chr+','+str(chr_pos+1)]
    multi_entry = entry.split('/')
    snp_list = []
    sample_list = []
    AF_list = []
    rsID_list = []
    snp_info_list = []
    for entry in multi_entry:
        split_entry = entry.split(';')
        sample_list.append(split_entry[0].split(','))
        snp_list.append(split_entry[1].split(',')[1])
        rsID_list.append(split_entry[2])
        AF_list.append(split_entry[3])
        snp_info_list.append(
            chr+'_'+str(chr_pos)+'_'+split_entry[1].split(',')[0]+'_'+split_entry[1].split(',')[1])
    return snp_list, sample_list, rsID_list, AF_list, snp_info_list


current_guide_chr_pos_direction = ''
listCluster = []
inTarget.readline()
for line in inTarget:
    split = line.strip().split('\t')
    guide_no_bulge = split[1].replace('-', '')
    # if (guide_no_bulge + split[3] + split[5] + split[6]) == current_guide_chr_pos_direction:
    realTarget = split[2]
    replaceTarget = split[2].replace('-', '')
    refSeq = genomeStr[int(split[4]):int(split[4])+len(replaceTarget)]
    # print(refSeq)
    # print(replaceTarget)
    replaceTargetsDict = dict()

    if split[6] == '-':
        replaceTarget = reverse_complement_table(replaceTarget)

    print('ref', refSeq)
    print('var original', replaceTarget)

    totalDict = dict()
    for pos_c, c in enumerate(replaceTarget):
        if c in iupac_code:
            snpToReplace, sampleSet, rsID, AF_var, snpInfo = retrieveFromDict(
                pos_c+int(split[4]))
            for i, elem in enumerate(snpToReplace):
                listReplaceTarget = list(refSeq)
                listReplaceTarget[pos_c] = elem
                totalDict[(pos_c, elem)] = dict()
                listInfo = [[rsID[i], AF_var[i], snpInfo[i]]]
                totalDict[(pos_c, elem)][0] = [listReplaceTarget,
                                               sampleSet[i], listInfo]

    # the time of the universe
    for key in totalDict:  # for each snp in target (fixpoint)
        count = 0  # to create levels of tuple (couples,triplets,...)
        # for each other snp in target (> fixpoint)
        for newkey in totalDict:
            if newkey[1] > key[1]:
                resultSet = totalDict[key][count][1].intersection(
                    totalDict[newkey][0][1])  # extract intersection of sample to generate possible multisnp target
                if len(resultSet) > 0:  # if set is not null
                    # go forward one level (couples,triplets,...)
                    count += 1
                    # add new snp to preceding target seq with snp
                    replaceTarget1 = totalDict[newkey][0][0]
                    replaceTarget2 = totalDict[key][count-1][0]
                    replaceTarget2[newkey] = replaceTarget1[newkey]
                    listInfo2 = totalDict[key][count-1][2].copy()
                    listInfo2.append(totalDict[newkey][0][2])
                    # add to new level the modified seq and set of samples and info of snp
                    totalDict[key][count] = [
                        replaceTarget2, resultSet, listInfo2]
                    # remove the new generated sample set from all lower levels
                    totalDict[key][count-1][1] = totalDict[key][count - 1][1] - \
                        totalDict[key][count][1]
                    totalDict[newkey][0][1] = totalDict[newkey][0][1] - \
                        totalDict[key][count][1]

    for key in totalDict:
        for level in totalDict[key]:
            if len(totalDict[key][level][1]) > 0:
                print(totalDict[key][level])
