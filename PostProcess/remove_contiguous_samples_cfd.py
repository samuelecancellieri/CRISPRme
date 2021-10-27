#!/usr/bin/env python
"""
Created on Fri Aug 28 15:58:04 2020

@author: francesco
"""
import sys
import time

'''
def get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info):
    list_ref = []
    dict_var = dict()
    for ele in cluster:
        if ele[snp_info] == 'n':
            list_ref.append(ele)
        else:
            if ele[snp_info] in dict_var.keys():
                dict_var[ele[snp_info]].append(ele)
            else:
                dict_var[ele[snp_info]] = [ele]
    
    list_ref.sort(key = lambda x : x[total])
    if len(list_ref) > 1:
        best_ref = list_ref[0]
        for ele_ref in list_ref[1:]:
            if float(ele_ref[cfd]) > float(best_ref[cfd]):
                fileOut_disc.write("\t".join(best_ref))
                best_ref = ele_ref
            else:
                fileOut_disc.write("\t".join(ele_ref))
        fileOut.write("\t".join(best_ref))
    elif len(list_ref) == 1:
        fileOut.write("\t".join(list_ref[0]))
    
    for key in dict_var.keys():
        list_var = dict_var[key]
        list_var.sort(key = lambda x : x[total])
        best_var = list_var[0]
        if len(list_var) > 1:
            for ele_var in list_var[1:]:
                if float(ele_var[cfd]) > float(best_var[cfd]):
                    fileOut_disc.write("\t".join(best_var))
                    best_var = ele_var
                else:
                    fileOut_disc.write("\t".join(ele_var))
            fileOut.write("\t".join(best_var))
        else:
            fileOut.write("\t".join(best_var))
'''


def get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info):
    final_list = []
    list_ref = []
    dict_var = dict()
    for ele in cluster:
        if ele[snp_info] == 'n':
            list_ref.append(ele)
        else:
            # merge samples of identical targets (coming from different VCF datasets)
            if (ele[pos], ele[snp_info]) in dict_var.keys():
                dict_var[(ele[pos], ele[snp_info])][0][true_guide - 2] = dict_var[(ele[pos], ele[snp_info])
                                                                                  ][0][true_guide - 2] + "," + ele[true_guide - 2]  # true_guide - 2 points to samples column
                dict_var[(ele[pos], ele[snp_info])][0][snp_info - 2] = dict_var[(ele[pos], ele[snp_info])
                                                                                ][0][snp_info - 2] + "," + ele[snp_info - 2]  # snp_info - 2 points to rsID column
                dict_var[(ele[pos], ele[snp_info])][0][snp_info - 1] = dict_var[(ele[pos], ele[snp_info])
                                                                                ][0][snp_info - 1] + "," + ele[snp_info - 1]  # ttuesnp_info_guide - 2 points to AF column
            else:
                dict_var[(ele[pos], ele[snp_info])] = [ele]

    list_ref.sort(key=lambda x: int(x[total]))
    var_only = False
    best_ref = ''
    if len(list_ref) > 1:
        best_ref = list_ref[0]
        for ele_ref in list_ref[1:]:
            if float(ele_ref[cfd]) > float(best_ref[cfd]):
                best_ref = ele_ref
        final_list.append(best_ref)
    elif len(list_ref) == 1:
        best_ref = list_ref[0]
        final_list.append(list_ref[0])
    else:
        var_only = True

    for key in dict_var.keys():
        list_var = dict_var[key]
        list_var.sort(key=lambda x: int(x[total]))
        best_var = list_var[0]
        if len(list_var) > 1:
            for ele_var in list_var[1:]:
                if float(ele_var[cfd]) > float(best_var[cfd]):
                    best_var = ele_var
            if var_only:
                best_var[12] = 'y'
            if var_only == False:
                # use best aligned ref seq in bestvar target
                best_var[3] = best_ref[3]
                # use best aligned ref seq score in bestvar target
                best_var[cfd+1] = best_ref[cfd+1]
            final_list.append(best_var)
        else:
            if var_only:
                best_var[12] = 'y'
            if var_only == False:
                # use best aligned ref seq in bestvar target
                best_var[3] = best_ref[3]
                # use best aligned ref seq score in bestvar target
                best_var[cfd+1] = best_ref[cfd+1]
            final_list.append(best_var)

    temp_final_list = list()
    for target in final_list:
        # remove duplicates into snp info col
        target[snp_info] = ','.join(set(target[snp_info].split(',')))
        # remove duplicate into rsID col
        target[snp_info-2] = ','.join(set(target[snp_info-2].split(',')))
        # remove duplicate into AF col
        target[snp_info-1] = ','.join(set(target[snp_info-1].split(',')))
        # remove duplicate into samples col
        target[true_guide-2] = ','.join(set(target[true_guide-2].split(',')))
        # append to temp list
        temp_final_list.append(target)

    # final list with polished targets (no duplicates in snp data)
    final_list = temp_final_list
    n_ele = len(final_list)
    # sort by total in ascending order
    # final_list.sort(key=lambda x: int(x[total]))
    if n_ele > 1:
        if str(final_list[0][cfd]) != '-1.0' and sort_order == 'score':
            # sort by score in descending order
            final_list.sort(key=lambda x: float(x[cfd]), reverse=True)
            # sort ascending total and then descending score
            # sorted_final_list = sorted(
            #     sorted(final_list, key=lambda x: int(x[total]), reverse=False), key=lambda x: float(x[cfd]), reverse=True)
            # final_list = sorted_final_list
            # select ref as besttarget to save
            if final_list[0][cfd] == final_list[1][cfd] and final_list[1][cfd-2] == 'n':
                final_list[1][cfd-1] = str(n_ele-1)  # bestSCORE
                fileOut.write("\t".join(final_list[1]))
                bestTarget = final_list.pop(1)
            else:  # select var as besttarget to save
                final_list[0][cfd-1] = str(n_ele-1)  # best SCORE
                fileOut.write("\t".join(final_list[0]))
                bestTarget = final_list.pop(0)
        else:
            # sort for total (mm+bul) in target
            final_list.sort(key=lambda x: int(x[total]))
            final_list[0][cfd-1] = str(n_ele-1)  # bestMM_BUL when not scored
            fileOut.write("\t".join(final_list[0]))
            bestTarget = final_list.pop(0)
        for ele in final_list:
            final_list[0][cfd-1] = str(n_ele-1)  # discarded targets in cluster
            fileOut_disc.write(("\t".join(ele)))
    elif n_ele == 1:
        fileOut.write("\t".join(final_list[0]))


tau = int(sys.argv[3])  # range in bp to merge targets
chrom = int(sys.argv[4])-1  # chromoso
pos = int(sys.argv[5])-1  # position of target
total = int(sys.argv[6])-1  # mm+bul value
true_guide = int(sys.argv[7])-1  # real guide used in the search
snp_info = int(sys.argv[8])-1  # snp_info (ref_alt_allele)
cfd = int(sys.argv[9])-1  # CFD score
sort_order = str(sys.argv[10])
# -1 is to get the correct "python enumeration" from the bash script

start = time.time()
with open(sys.argv[1], 'r') as fileIn:
    header = fileIn.readline()
    with open(sys.argv[2], 'w') as fileOut:
        with open(sys.argv[2]+'.discarded_samples', 'w') as fileOut_disc:
            fileOut.write(header)
            fileOut_disc.write(header)
            prev_pos = -(tau+1)
            best_row = ""
            prev_guide = ""
            prev_chr = ""
            prev_snp = ""
            cluster = []
            for line in fileIn:
                splitted = line.split("\t")
                if prev_guide != splitted[true_guide] or prev_chr != splitted[chrom] or int(splitted[pos]) - prev_pos > tau:
                    get_best_targets(cluster, fileOut,
                                     fileOut_disc, cfd, snp_info)
                    cluster = [splitted]
                else:
                    cluster.append(splitted)
                prev_guide = splitted[true_guide]
                prev_pos = int(splitted[pos])
                prev_chr = splitted[chrom]
                prev_snp = splitted[snp_info]

            get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info)

print("Mergin done in: "+str(time.time()-start))
