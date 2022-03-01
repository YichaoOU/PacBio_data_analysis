#!/usr/bin/env python3

import sys


fileName = open("temp_SVIM.vcf", 'r')
Sniffles_file = open("temp_sniffles.vcf", 'r')
outFile = open('temp_no_repeated_sniffles_svim.vcf', 'w+')

SUP_READS_SVIM = int(sys.argv[2])


def make_reads_unify(read_name_vec):
    vec_names = {}
    for i, vec in enumerate(read_name_vec):
        # name = vec[vec.find('=') + 1: vec.find(';')]
        if vec not in vec_names.keys():
            vec_names[vec] = 1
        else:
            vec_names[vec] += 1

    return vec_names


def reconstruct_reads(vec_names, line7):
    read_names = ",".join(list(vec_names.keys()))
    if line7.find(";SUPTYPE") != -1:
        line7 = line7.replace(line7[line7.find("READS="):line7.find(";SUPTYPE")], "READS="+read_names)
    else:
        line7 = line7.replace(line7[line7.find("READS="):], "READS="+read_names)

    return line7


for line in Sniffles_file:
    outFile.write(line)

for line in fileName:
    if line[0] != '#':
        line_vec = line.split('\t')
        if line_vec[7].find("RNAMES=") != -1:
            line_vec[7] = line_vec[7].replace(line_vec[7][line_vec[7].find("RNAMES"):line_vec[7].find("RNAMES")+7]
                                              , "READS=")
            read_name_vec = line_vec[7].split("READS=")
            read_name_vec[1] = read_name_vec[-1].replace(read_name_vec[-1][read_name_vec[-1].find(";SUPTYPE"):], "")
            just_reads = read_name_vec[1].split(",")
        else:
            read_name_vec = line_vec[7].split("READS=")
            just_reads = read_name_vec[1].split(",")
        dict_name = make_reads_unify(just_reads)

        for val in dict_name.values():
            if val > 1:
                if len(dict_name) >= SUP_READS_SVIM:
                    line_vec[7] = reconstruct_reads(dict_name, line_vec[7])
                    break
        if len(dict_name) >= SUP_READS_SVIM:
            temp = '\t'.join(line_vec)
            outFile.write(temp)

fileName.close()
outFile.close()
