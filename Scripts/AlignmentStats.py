#this code opens iterates through aligned files and creates csv containing number of taxa, GC content, length and gap content
#averages for all alignments and standard deviations were then produced in RStudio


from Bio import SeqIO
import csv
import os
import statistics

from scipy import mean
from scipy.stats import sem, t


directory = "."

with open('../alignments_stats29.csv', mode = 'w') as align_file:
	stat_writer = csv.writer(align_file, delimiter=',', quotechar = '"', quoting=csv.QUOTE_MINIMAL)
	stat_writer.writerow(["ogname", "numtaxa", "avlength", "stdlength", "avGC", "stdGC", "avgap"])
	for file in os.listdir(directory):
		record = list(SeqIO.parse(open(file), 'fasta'))
		GC_total = 0
		length_total = 0
		gap_total = 0
		all_dna_lengths = []
		all_GC_content = []
		all_gap_content = []
		for sequence in record:
			dna_count = sequence.seq.count("a") + sequence.seq.count("t") + 	
				sequence.seq.count("g") + sequence.seq.count("c")
			dna_length = dna_count + sequence.seq.count("-")
			GC_num = (sequence.seq.count("g") + sequence.seq.count("c"))
			GC_ratio = GC_num / float(dna_count)
			length_total = length_total + dna_length
			GC_total = GC_total + GC_ratio
			gap_ratio = sequence.seq.count("-") / float(dna_length)
			all_dna_lengths.extend([dna_length])
			all_GC_content.extend([GC_total])
			all_gap_content.extend([gap_ratio])
		if length_total != 0:
			taxa = len(record)
			av_length = length_total/float(taxa)
			av_GC = GC_total/float(taxa)
			av_gap = gap_ratio/float(taxa)
			stdlength = statistics.stdev(all_dna_lengths)
			stdGC = statistics.stdev(all_GC_content)
		base = os.path.basename(file)
		filename = base.rsplit('_'[0])
		ogname = "{}_{}".format(filename[0], filename[1])
		print(ogname, stdlength, stdGC)
		stat_writer.writerow([ogname, taxa, av_length, stdlength, av_GC, stdGC, av_gap])
