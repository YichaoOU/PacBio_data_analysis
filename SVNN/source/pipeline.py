import os
import subprocess
import sys


def minimap2(reads, ref, thread=4):
	command = "../bin/minimap2 -t " + str(thread) + " --MD -ax map-hifi " + ref + " " + reads + " > " + "temp_mini.sam"
	print (command)
	subprocess.run(command, shell=True)


def find_sr():
	subprocess.run("../bin/find_SR")


def find_informative_reads():
	subprocess.run("python3 Neural_Net.py", shell=True)


def find_in_fastq(read_dir):
	command = "../bin/find_SA_in_FastQ " + read_dir
	subprocess.run(command, shell=True)


def NGMLR(ref, thread=4):
	command = "../bin/ngmlr --bam-fix -t " + str(thread) + " -r " + ref + " -q  temp_reads_to_ngmlr.fastq -o temp_ngmlr_just_split.sam"
	subprocess.run(command, shell=True)


def combine_NGMLR_minimap():
	command = "../bin/combined temp_mini.sam"
	subprocess.run(command, shell=True)


def samtools(thread=4):
	command = "../bin/samtools sort -o temp_mini_ngmlr.bam -@ " + str(thread) + " temp_mini_ngmlr.sam"
	subprocess.run(command, shell=True)


def sniffles(supporting_reads, thread=4):
	command = "../bin/sniffles --genotype -n -1 -t " + str(thread) + " -s " + str(supporting_reads) + \
			" -m temp_mini_ngmlr.bam -v temp_sniffles.vcf"
	subprocess.run(command, shell=True)


def svim(ref, supporting_reads):
	command = "../bin/svim alignment --skip_genotyping --read_names temp_SVIM temp_mini_ngmlr.bam " + ref
	subprocess.run(command, shell=True)
	command = "python3 filter_SV.py " + str(supporting_reads)
	subprocess.run(command, shell=True)


def clean_redundacny(sniffles_sup, svim_sup):
	command = "python3 unify_Supporting_Reads.py " + str(sniffles_sup) + " " + str(svim_sup)
	subprocess.run(command, shell=True)
	subprocess.run("python3 unify_Sniffles_SVIM.py", shell=True)


def make_standard(parameters):
	"""
	Makes arguments in a standard format so for example '-' and '–' will be the same
	"""
	for i in range(len(parameters)):
		if '–' in parameters[i]:
			parameters[i] = '-' + parameters[i][1:]
		if ',' in parameters[i]:
			parameters[i] = parameters[i][:-1]
	return parameters


def main():
	print('Argument List:', sys.argv)
	parameters = make_standard(sys.argv[1:])

	try:
		fastQ1 = parameters[parameters.index('-q')+1]
	except ValueError:
		raise ValueError("Read file not given")

	try:
		reference = parameters[parameters.index('-r')+1]
	except ValueError:
		raise ValueError("Reference not given")

	try:
		support_reads = parameters[parameters.index('-s1')+1]
	except ValueError:
		raise ValueError("Number of supporting reads for Sniffles not given")
	try:
		svim_supp_reads = parameters[parameters.index('-s2')+1]
	except ValueError:
		raise ValueError("Number of supporting reads for SVIM not given")

	print ("running minimap2")
	minimap2(fastQ1, reference)
	print ("running find_sr")
	find_sr()
	print ("running find_informative_reads")
	find_informative_reads()
	print ("running find_in_fastq")
	find_in_fastq(fastQ1)
	NGMLR(reference)
	combine_NGMLR_minimap()
	samtools()
	sniffles(support_reads)
	svim(reference, svim_supp_reads)
	clean_redundacny(support_reads, svim_supp_reads)

	subprocess.run("rm -rf temp_SVIM", shell=True)
	subprocess.run("rm temp_*", shell=True)


if __name__ == '__main__':
	main()

