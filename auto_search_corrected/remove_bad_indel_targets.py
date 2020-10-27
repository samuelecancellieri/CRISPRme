#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 13:55:37 2020

@author: francesco
"""
import sys
import os

with open(sys.argv[1], 'r') as fileIn:
    with open(sys.argv[1]+'.tmp', 'w') as fileOut:
        header = fileIn.readline()
        fileOut.write(header)
        for line in fileIn:
            splitted = line.split('\t')
            if splitted[18] != 'n':
                indel_info = splitted[18].split('_')
                start_indel = int(indel_info[1])
                if len(indel_info[2]) < len(indel_info[3]):
                    range_indel = range(start_indel, start_indel+len(indel_info[3]))
                else:
                    range_indel = range(start_indel, start_indel+len(indel_info[2]))
                range_target = range(int(splitted[5]), int(splitted[5])+len(splitted[1]))
                if set(range_target).intersection(range_indel):
                    fileOut.write(line)
            else:
                fileOut.write(line)

os.system('mv '+sys.argv[1]+'.tmp '+sys.argv[1])