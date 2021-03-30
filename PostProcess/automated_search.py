#!/usr/bin/env python

import sys
import os

script_path = os.path.dirname(os.path.abspath( __file__ ))

input_args = sys.argv

variant = True

if "--help" in input_args:
    print("This is the automated search process that goes from raw input up to the post-analysis of results.")
    print("These are the flags that must be used in order to run this function:")
    print("\t--genome , used to specify the reference genome folder")
    print("\t--vcf , used to specify the VCF folder [OPTIONAL!]")
    print("\t--guide , used to specify the file that contains guides used for the search")
    print("\t--pam , used to specify the file that contains the pam")
    print("\t--annotation , used to specify the file that contains some annotations of the reference genome")
    print("\t--samplesID , used to specify the file that contains the information about samples present in VCF files [OPTIONAL!]")
    print("\t--bMax , used to specify the number of bulges for the indexing of the genome(s)")
    print("\t--mm , used to specify the number of mismatches permitted in the search phase")
    print("\t--bDNA , used to specify the number of DNA bulges permitted in the search phase")
    print("\t--bRNA , used to specify the number of RNA bulges permitted in the search phase")
    print("\t--merge, used to specify the threshold used to merge close targets")
    print("\t--output , used to specify the output folder for the results")
    exit(0)

if "--genome" not in input_args:
    print("--genome must be contained in the input")
    exit(1)
else:
    try:
        genomedir = os.path.abspath(input_args[input_args.index("--genome")+1])
    except IndexError:
        print("Please input some parameter for flag --genome")
        exit(1)
    if not os.path.isdir(genomedir):
        print("The folder specified for --genome does not exist")
        exit(1)
    
if "--vcf" not in input_args: 
    variant = False
else:
    try:
        vcfdir = os.path.abspath(input_args[input_args.index("--vcf")+1])
    except IndexError:
        print("Please input some parameter for flag --vcf")
        exit(1)
    if not os.path.isdir(vcfdir):
        print("The folder specified for --vcf does not exist")
        exit(1)
    
if "--guide" not in input_args:
    print("--guide must be contained in the input")
    exit(1)
else:
    try:
        guidefile = os.path.abspath(input_args[input_args.index("--guide")+1])
    except IndexError:
        print("Please input some parameter for flag --guide")
        exit(1)
    if not os.path.isfile(guidefile):
        print("The folder specified for --guide does not exist")
        exit(1)
    
if "--pam" not in input_args:
    print("--pam must be contained in the input")
    exit(1)
else:
    try:
        pamfile = os.path.abspath(input_args[input_args.index("--pam")+1])
    except IndexError:
        print("Please input some parameter for flag --pam")
        exit(1)
    if not os.path.isfile(pamfile):
        print("The folder specified for --pam does not exist")
        exit(1)
    
if "--annotation" not in input_args:
    print("--annotation must be contained in the input")
    exit(1)
else:
    try:
        annotationfile = os.path.abspath(input_args[input_args.index("--annotation")+1])
    except IndexError:
        print("Please input some parameter for flag --annotation")
        exit(1)
    if not os.path.isfile(annotationfile):
        print("The folder specified for --annotation does not exist")
        exit(1)
    
if variant and "--samplesID" not in input_args:
    print("--samplesID must be contained in the input")
    exit(1)
elif not variant and "--samplesID" in input_args:
    print("--samplesID was in the input but no VCF directory was specified")
    exit(1)
elif "--samplesID" in input_args:
    try:
        samplefile = os.path.abspath(input_args[input_args.index("--samplesID")+1])
    except IndexError:
        print("Please input some parameter for flag --samplesID")
        exit(1)
    if not os.path.isfile(samplefile):
        print("The folder specified for --samplesID does not exist")
        exit(1)
    
if "--bMax" not in input_args:
    print("--bMax must be contained in the input")
    exit(1)
else:
    try:
        bMax = input_args[input_args.index("--bMax")+1]
    except IndexError:
        print("Please input some parameter for flag --bMax")
        exit(1)
    try:
        bMax = int(bMax)
    except:
        print("Please input a number for flag bMax")
        exit(1)
    if bMax < 0 or bMax > 2:
        print("The range for bMax is from 0 to 2")
        exit(1)

if "--mm" not in input_args:
    print("--mm must be contained in the input")
    exit(1)
else:
    try:
        mm = input_args[input_args.index("--mm")+1]
    except IndexError:
        print("Please input some parameter for flag --mm")
        exit(1)
    try:
        mm = int(mm)
    except:
        print("Please input a number for flag mm")
        exit(1)
    
if "--bDNA" not in input_args:
    print("--bDNA must be contained in the input")
    exit(1)
else:
    try:
        bDNA = input_args[input_args.index("--bDNA")+1]
    except IndexError:
        print("Please input some parameter for flag --bDNA")
        exit(1)
    try:
        bDNA = int(bDNA)
    except:
        print("Please input a number for flag bDNA")
        exit(1)
    if bDNA > bMax:
        print("The number of bDNA must be equal or less than bMax")
        exit(1)
    elif bDNA < 0 or bDNA > 2:
        print("The range for bDNA is from 0 to", bMax)
        exit(1)
    
if "--bRNA" not in input_args:
    print("--bRNA must be contained in the input")
    exit(1)
else:
    try:
        bRNA = input_args[input_args.index("--bRNA")+1]
    except IndexError:
        print("Please input some parameter for flag --bRNA")
        exit(1)
    try:
        bRNA = int(bRNA)
    except:
        print("Please input a number for flag bRNA")
        exit(1)
    if bRNA > bMax:
        print("The number of bRNA must be equal or less than bMax")
        exit(1)
    elif bRNA < 0 or bRNA > 2:
        print("The range for bRNA is from 0 to", bMax)
        exit(1)

if "--merge" not in input_args:
    print("--merge must be contained in the input")
    exit(1)
else:
    try:
        merge_t = input_args[input_args.index("--merge")+1]
    except IndexError:
        print("Please input some parameter for flag --merge")
        exit(1)
    try:
        merge_t = int(merge_t)
    except:
        print("Please input a number for flag merge")
        exit(1)
    if merge_t < 0:
        print("Please specify a positive number for --merge")
        exit(1)
    
if "--output" not in input_args:
    print("--output must be contained in the input")
    exit(1)
else:
    try:
        outputfolder = os.path.abspath(input_args[input_args.index("--output")+1])
    except IndexError:
        print("Please input some parameter for flag --output")
        exit(1)
    if not os.path.isdir(outputfolder):
        print("The folder specified for --output does not exist")
        exit(1)
    
os.chdir(script_path)
if variant:
    os.system("./automated_search_good_parallel_v2.sh "+genomedir+" "+vcfdir+" "+guidefile+" "+pamfile+" "+annotationfile+" "+samplefile+" "+str(bMax)+" "+str(mm)+" "+str(bDNA)+" "+str(bRNA)+" "+str(merge_t)+" "+outputfolder+" "+script_path+" "+str(len(os.sched_getaffinity(0))-2))
else:
    os.system("./automated_search_good_parallel_v2.sh "+genomedir+" _ "+guidefile+" "+pamfile+" "+annotationfile+" _ "+str(bMax)+" "+str(mm)+" "+str(bDNA)+" "+str(bRNA)+" "+str(merge_t)+" "+outputfolder+" "+script_path+" "+str(len(os.sched_getaffinity(0))-2))
    
    