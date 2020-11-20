"""
Created on Fri Aug 28 15:58:04 2020

@author: francesco
"""
import sys
import time

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
    
    

tau = int(sys.argv[3])
chrom = int(sys.argv[4])-1 
pos = int(sys.argv[5])-1 
snp_info = int(sys.argv[6])-1 
cfd = int(sys.argv[7])-1
# -1 is to get the correct "python enumeration" from the bash script

start = time.time()
with open(sys.argv[1], 'r') as fileIn:
    header = fileIn.readline()
    with open(sys.argv[2], 'w') as fileOut:
        with open(sys.argv[2]+'.discarded_samples', 'w') as fileOut_disc:
            fileOut.write(header)
            prev_pos = -(tau+1)
            best_row = ""
            prev_chr = ""
            prev_snp = ""
            cluster = []
            for line in fileIn:
                splitted = line.split("\t")
                if prev_chr != splitted[chrom] or int(splitted[pos]) - prev_pos > tau:
                    get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info)
                    cluster = [splitted]
                else:
                    cluster.append(splitted)
                prev_pos = int(splitted[pos])
                prev_chr = splitted[chrom]
                prev_snp = splitted[snp_info]
            
            get_best_targets(cluster, fileOut, fileOut_disc, cfd, snp_info)

print("Mergin done in: "+str(time.time()-start))
