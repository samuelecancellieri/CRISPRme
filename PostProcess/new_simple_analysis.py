#!/usr/bin
from os import replace
import sys
import json

inFasta = open(sys.argv[1], 'r')  # lettura fasta del chr
inFasta.readline()  # lettura fasta del chr
genomeStr = inFasta.readlines()  # lettura fasta del chr
genomeStr = ''.join(genomeStr).upper()
genomeStr = genomeStr.replace('\n', '')
inDict = sys.argv[2]
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


def retrieveFromDict():
    # bla
    print('ciao')


current_guide_chr_pos_direction = ''
listCluster = []
for line in inTarget:
    split = line.strip().split('\t')
    guide_no_bulge = split[1].replace('-', '')
    if (guide_no_bulge + split[3] + split[5] + split[6]) == current_guide_chr_pos_direction:
        realTarget = split[2]
        replaceTarget = split[2].replace('-', '')
        refSeq = genomeStr[split[4]:split[4]+len(guide_no_bulge)]
        replaceTargetsDict = dict()

        totalDict = dict()
        for pos_c, c in enumerate(replaceTarget):
            if c in iupac_code:
                snpToReplace, sampleSet, rsID, AF_var, snpInfo = retrieveFromDict(
                    pos_c+int(split[4]))
                listReplaceTarget = list(refSeq)
                listReplaceTarget[pos_c] = snpToReplace
                totalDict[pos_c] = dict()
                listInfo = [[rsID, AF_var, snpInfo]]
                totalDict[pos_c][0] = [listReplaceTarget,
                                       sampleSet, listInfo]

        # the time of the universe
        for key in totalDict:  # for each snp in target (fixpoint)
            count = 0  # to create levels of tuple (couples,triplets,...)
            # for each other snp in target (> fixpoint)
            for newkey in totalDict:
                if newkey > key:
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
