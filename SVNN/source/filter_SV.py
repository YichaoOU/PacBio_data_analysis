#!/usr/bin/env python3

import sys


fileName = open("temp_SVIM/variants.vcf", 'r')
outFile = open('temp_SVIM.vcf', 'w+')
min_count = int(sys.argv[1])


for line in fileName:
    if line[0] != '#':
        line_vec = line.split('\t')

        if line_vec[7].count(",") + 1 >= min_count:
            outFile.write(line)
    else:
        outFile.write(line)

fileName.close()
outFile.close()        
