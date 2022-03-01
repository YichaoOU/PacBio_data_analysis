#!/usr/bin/env python3

inFile = open("temp_no_repeated_sniffles_svim.vcf", 'r')
outFile = open("final_results.vcf", "w+")

listOfLines = inFile.readlines()

MAX_DIFF = 50
chr21_dict = {}


def extract_info(line_vec):
    chr_num = line_vec[0]
    start = line_vec[1]
    type_str_idx = line_vec[7].find('SVTYPE=')
    type_end_idx = line_vec[7].find(';', type_str_idx)
    type_sv = line_vec[7][type_str_idx + 7:type_end_idx]
    if type_sv == 'DUP:TANDEM':
        type_sv = 'DUP'

    if type_sv == 'BND':
        end_str_idx = line_vec[4].find(':')
        end_end_idx = max(line_vec[4].find('[', end_str_idx), line_vec[4].find(']', end_str_idx))
        end = line_vec[4][end_str_idx + 1:end_end_idx]
    else:
        end_str_idx = line_vec[7].find('END=')
        end_end_idx = line_vec[7].find(';', end_str_idx)
        end = line_vec[7][end_str_idx + 4:end_end_idx]

    return chr_num, type_sv, start, end


line_mem = []
for i, line in enumerate(listOfLines):
    if i not in line_mem:
        outFile.write(line)
        if line[0] != '#':
            line_vec = line.split('\t')
            chr_num1, type_sv1, start1, end1 = extract_info(line_vec)
            for j, line_2nd in enumerate(listOfLines[i+1:]):
                line_2nd_vec = line_2nd.split('\t')
                chr_num2, type_sv2, start2, end2 = extract_info(line_2nd_vec)
                if (abs(int(start1) - int(start2)) < MAX_DIFF) and (abs(int(end1) - int(end2)) < MAX_DIFF):
                    if chr_num1 == chr_num2:
                        line_mem.append(i+j+1)
