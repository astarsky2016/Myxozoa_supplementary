from Bio import SeqIO
from collections import defaultdict
import pickle, sys, csv, plyvel
from ete3 import NCBITaxa
tax_data = {}
calc = []
ncbi = NCBITaxa()
def getSeq():
	scan_results = defaultdict(list)
	with open("SILVA/ssrna_scan_position_assembly.tab") as s_tab:
		for line in s_tab:
			data = line.split()
			ctg = data[0]
			ids = data[1]
			start = int(data[2])
			end = int(data[3])
			med = int((start + end)/2)
			e_val = float(data[4])
			if (e_val < 1e-10) and (abs(end-start) > 400):
				scan_results[ctg].append([ids, start, end, e_val])
	contigs = {}
	for rec in SeqIO.parse("flye_single_run1/assembly.fasta", "fasta"):
		contigs[rec.id] = rec
	print ("loaded")
	ribosomal = []
	for ctg in scan_results.keys(): 
		sorted_res = sorted(scan_results[ctg], key=lambda tup: tup[-1])	
		limits = [0]
		pos = 0
		for candidate in sorted_res:
			start = candidate[1]
			end = candidate[2]
			if start >= limits[-1]:
				pos += 1
				limits.append(end)
				candidate_ribosomal = contigs[ctg][start:end]
				candidate_ribosomal.id = "{}_{}_18SrRNA".format(ctg, pos)
				ribosomal.append(candidate_ribosomal)
				print(ctg, pos)
	SeqIO.write(ribosomal, open("18S/rRNA_candidates.fasta", "w"), "fasta")
	print (len(ribosomal))
def screenCsv():
	db = plyvel.DB('./acc2tax_levDB_1')
	results = defaultdict(lambda: defaultdict(int))
	with open("cnidaria.csv") as csv_file:
		for row in csv.reader(csv_file):
			contig = "_".join(row[0].split("_")[0:2])
			acc = row[1]
			percent = float(row[2])
			aln_len = int(row[3])
			e_val = float(row[-2])
			if ((percent > 90.0) and (e_val < 1e-10) and (aln_len > 400)):
				t = db.get(bytes(acc, 'utf8'))
				if t:
					results[contig][int(t)] += 1
					print ("dodajem")
	db.close()
	print ("parsed")
	f_out = open("krona_taxa.txt", "w")
	for ctg in results.keys():
		for animal in results[ctg].keys():
			lineage = ncbi.get_lineage(animal)
			taxid2name = ncbi.get_taxid_translator(lineage)
			named = "\t".join([taxid2name[l] for l in lineage])
			f_out.write("{}\t".format(results[ctg][animal]) + named + "\n")
	f_out.close()
def processKrona():
	for rec in SeqIO.parse("SILVA/SILVA_138.1_SSURef_NR99_tax_silva.fasta", "fasta"):
		ids = rec.id
		taxonomy = rec.description.split()[1].strip()
		if "Eukaryota" in taxonomy:
			calc.append(len(rec.seq))
		tax_data[ids] = taxonomy
	print (max(calc), min(calc), sum(calc)/len(calc))
	print ("taxonomy mapped")
	scan_results = defaultdict(list)
	with open("SILVA/ssrna_scan_position.tab") as s_tab:
		for line in s_tab:
			data = line.split()
			ctg = data[0]
			ids = data[1]
			start = int(data[2])
			end = int(data[3])
			med = int((start + end)/2)
			e_val = float(data[4])
			if e_val < 1e-10:
				scan_results[ctg].append([ids, abs(end-start), med, e_val])
	print ("data processed")
	avg_sum = []
	for ctg in scan_results.keys(): 
		sorted_res = sorted(scan_results[ctg], key=lambda tup: tup[-1])
		avg_sum.append(sorted_res[0][1])
	print (max(avg_sum))
	avg_sum = int(sum(avg_sum)/len(avg_sum))
	print (avg_sum)
	#tolerance = 10
	cnt = 0
	krona_input = defaultdict(int)
	for ctg in scan_results.keys(): 
		data = scan_results[ctg]
		positional = [0]
		for match in data:
			if match[1] > 400:
				cnt += 1
				for p in positional:
					if abs(match[2] - p) > 2*avg_sum:
						positional.append(match[2])
						krona_input[tax_data[match[0]]] += 1
	pickle.dump(dict(krona_input), open("krona_tax.p", "wb"))
	with open("krona_taxa.txt", "w") as f_in:
		for t in krona_input.keys():
			count = krona_input[t]
			f_in.write("{}".format(count) + "\t" + "\t".join(t.split(";")) + "\n")
if __name__ == '__main__':
	getSeq()





	
