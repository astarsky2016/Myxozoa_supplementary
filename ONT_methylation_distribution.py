from Bio import SeqIO
import sys, plyvel, pickle, re, itertools, glob, os
from collections import defaultdict
import numpy as np
from scipy.stats import kstest
def find_CpG(genome):
	CpGs = defaultdict(list)
	pattern = r"(CG|GC)"
	total = 0
	seq_len = 0
	chromosome_sizes = {}
	for rec in SeqIO.parse(genome, "fasta"):
		chromosome_sizes[rec.id] = len(rec.seq)
		matches = re.finditer(pattern, str(rec.seq).upper())
		# Print the positions of CpG sites
		for match in matches:
			CpGs[rec.id].append(match.start())
			total += 1
	print (total)
	for ctg in CpGs.keys():
		data = sorted(CpGs[ctg])
		# Expected cumulative distribution function (CDF) for a uniform distribution
		uniform_cdf = np.around((np.arange(1, len(data) + 1) / len(data)) * chromosome_sizes[ctg]).astype(int)
		# Perform the KS test
		ks_statistic, p_value = kstest(data, 'uniform', args=(0, len(data)))
		# Check the p-value
		alpha = 0.05  # Significance level
		if p_value < alpha:
			#print("The data does not follow a uniform distribution (reject null hypothesis).")
			pass
		else:
			print (ks_statistic, p_value)
			print("The data follows a uniform distribution (fail to reject null hypothesis).")
	return dict(CpGs)
def multiGenome():
	data = defaultdict(list)
	for gmark_in in glob.glob("/disk3/ALL_ASSEMBLIES/ORIGINAL_ASSEMBLY/*_metaeuk.gff"):
		assembly = gmark_in.replace("_metaeuk.gff", ".fna")
		name = os.path.basename(gmark_in).replace("_metaeuk.gff", "")
		gmark = genemark(gmark_in)
		stats = GC(gmark, assembly)
		data["Species"].append(name)
		data["genome GC"].append(stats["genome"][0])
		data["genome CpG O/E"].append(stats["genome"][1])
		data["CDS GC"].append(stats["CDS"][0])
		data["CDS CpG O/E"].append(stats["CDS"][1])
		data["Non-CDS GC"].append(stats["Non-CDS"][0])
		data["Non-CDS CpG O/E"].append(stats["Non-CDS"][1])		
	pickle.dump(dict(data), open("/disk3/ALL_ASSEMBLIES/GC_data.p", "wb"))

def genome(f):
	total = 0
	patterns = [r"({}{})".format(el[0], el[1]) for el in itertools.permutations('ACGT', 2)]
	counts = [0]*len(patterns)
	gc = 0
	for rec in SeqIO.parse(f, "fasta"):
		total += len(rec.seq)
		dna = str(rec.seq).upper()
		for base in dna:
			if base == "G" or base == "C":
				gc += 1
		for p, pattern in enumerate(patterns):
			matches = [match.start() for match in re.finditer(pattern, dna)]
			counts[p] += len(matches)*2
	print ("GC: " "{:.2f}%".format((gc/total)*100))
	print ([(patterns[pos], "{:.1f}".format((el/total)*100)) for pos, el in enumerate(counts)])
def find_sequential_parts(input_list):
    sequential_parts = []
    current_part = [input_list[0]]

    for i in range(1, len(input_list)):
        if input_list[i] == input_list[i - 1] + 1:
            current_part.append(input_list[i])
        else:
            sequential_parts.append(current_part)
            current_part = [input_list[i]]

    sequential_parts.append(current_part)  # Add the last part
    GC_count = 0
    for el in sequential_parts:
    	if len(el) == 1:
    		GC_count += 2
    	else:
    		GC_count += len(el) + 1
    return GC_count
def all_sample_genes():
	gene_content = []
	pattern = r"(CG)"
	samples = glob.glob("ORIGINAL_ASSEMBLY/*.codon.fas")
	for genes in samples:
		for rec in SeqIO.parse(genes, "fasta"):
			data = rec.id.split("|")
			acc = data[0]
			ctg = data[1]
			dna = str(rec.seq).upper()
			matches = [match.start() for match in re.finditer(pattern, dna)]
			gene_content.append((len(matches), len(dna), acc, ctg))
		print (genes)
	pickle.dump(gene_content, open("all_sample_gene_CpG_content.p", "wb"))
def gene_CpG_all():
	db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel/')
	pattern = r"(CG)"
	gene_content = []
	gc = 0
	gcat = 0
	cpg = 0
	acc2func = defaultdict(list)
	for gene in glob.glob("ORIGINAL_ASSEMBLY/*.codon.fas"):
		if "sample" in gene:
			basename = os.path.basename(gene).replace("_metaeuk.codon.fas", "")
			print (gene)
			for rec in SeqIO.parse(gene, "fasta"):
				dna = str(rec.seq).upper()
				count_g = 0
				count_c = 0
				for base in dna:
					if base == "G": count_g += 1
					if base == "C": count_c += 1
					if base == "G" or base == "C":
						gc += 1
					gcat += 1
				matches = [match.start() for match in re.finditer(pattern, dna)]
				cg_num = len(matches)
				try:
					CpG_OE = (cg_num*len(dna))/(count_c*count_g)
				except:
					CpG_OE = None
				rec_data = rec.id.split("|")
				acc = rec_data[0]
				ctg = rec_data[1]
				start = int(rec_data[6])
				end = int(rec_data[7])
				go = db.get(acc.encode('utf-8'))
				if go is not None:
					res = go.decode('utf-8').split("|")
					for r in res:
						desc = r.split("#")[1]
						desc_data = desc.split(":")
						category = desc_data[0]
						desc = desc_data[1]
						acc2func[acc].append(desc_data)
				gene_content.append((len(matches), len(dna), basename, acc, CpG_OE))
	pickle.dump(gene_content, open("gene_CpG_test_ONT.p", "wb"))
	pickle.dump(acc2func, open("test_genes_desc_ONT.p", "wb"))
def gene_CpG_test(genes):
	db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel/')
	pattern = r"(CG)"
	gene_content = []
	gc = 0
	gcat = 0
	cpg = 0
	acc2func = defaultdict(list)
	for gene in glob.glob("ORIGINAL_ASSEMBLY/*.codon.fas"):
		if "sample" in gene:
			basename = os.path.basename(gene).replace("_metaeuk.codon.fas", "")
			meth = gene.replace("_metaeuk.codon.fas", "_methylation_frequency.tsv")
			meth_pos = mapMeth(meth)
			for rec in SeqIO.parse(gene, "fasta"):
				dna = str(rec.seq).upper()
				count_g = 0
				count_c = 0
				for base in dna:
					if base == "G": count_g += 1
					if base == "C": count_c += 1
					if base == "G" or base == "C":
						gc += 1
					gcat += 1
				matches = [match.start() for match in re.finditer(pattern, dna)]
				cg_num = len(matches)
				try:
					CpG_OE = (cg_num*len(dna))/(count_c*count_g)
				except:
					CpG_OE = None
				rec_data = rec.id.split("|")
				acc = rec_data[0]
				ctg = rec_data[1]
				start = int(rec_data[6])
				end = int(rec_data[7])
				meth_sum = []
				go = db.get(acc.encode('utf-8'))
				if go is not None:
					res = go.decode('utf-8').split("|")
					for r in res:
						desc = r.split("#")[1]
						desc_data = desc.split(":")
						category = desc_data[0]
						desc = desc_data[1]
						acc2func[acc].append(desc_data)
				for coord in meth_pos.get(ctg, []):
					if coord[0] in range(start, end + 1) and coord[1] in range(start, end + 1):
						meth_sum.append(coord[2])
				if meth_sum:
					meth_cumulative = sum(meth_sum)/len(meth_sum)
				else:
					meth_cumulative = None
				gene_content.append((len(matches), len(dna), basename, rec.id, meth_cumulative, CpG_OE))
	pickle.dump(gene_content, open("gene_CpG_test.p", "wb"))
	#pickle.dump(acc2func, open("test_genes_desc.p", "wb"))
def gene_CpG(genes):
	#core_func = pickle.load(open("core_GO.p", "rb"))
	db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel/')
	#pattern = r"(CG|GC)"
	pattern = r"(CG)"
	gene_content = []
	core_genes = set()
	gc = 0
	gcat = 0
	cpg = 0
	acc2func = defaultdict(list)
	for rec in SeqIO.parse(genes, "fasta"):
		data = rec.id.split("|")
		acc = data[0]
		ctg = data[1]
		dna = str(rec.seq).upper()
		for base in dna:
			if base == "G" or base == "C":
				gc += 1
			gcat += 1
		matches = [match.start() for match in re.finditer(pattern, dna)]
		#cpg += find_sequential_parts(matches)
		cpg += 2*len(matches)
		gene_content.append((len(matches), len(dna), acc, ctg))
		go = db.get(acc.encode('utf-8'))
		if go is not None:
			res = go.decode('utf-8').split("|")
			for r in res:
				desc = r.split("#")[1]
				desc_data = desc.split(":")
				category = desc_data[0]
				desc = desc_data[1]
				core_genes.add(acc)
				acc2func[acc].append(desc_data)
				#if category == "P":
				#	core_desc = core_func[category]
				#	if desc in core_desc:
				#		core_genes.add(acc)
	pickle.dump(gene_content, open("gene_CpG_content.p", "wb"))
	print (len(core_genes))
	pickle.dump(core_genes, open("core_genes_acc.p", "wb"))
	pickle.dump(acc2func, open("core_genes_desc.p", "wb"))
	print ("CpG role:", gcat, gc/gcat, cpg/gcat)
	print(cpg)
	print ((cpg/gcat)/((gc/gcat)/2)**2)
def genemark(f):
	contig_genes = defaultdict(list)
	lns = []
	total = 0
	with open(f, "r") as f_in:
		for line in f_in:
			data = [el.strip() for el in line.split()]
			if data[2] == "gene":
				start = int(data[3])
				end = int(data[4])
				ctg = data[0]
				size = (end - start) + 1
				lns.append(size)
				if size >= 50: #300 za GeneMark, 50 za MetaEuk
					contig_genes[ctg].append((start, end))
					total += size
	print (min(lns), max(lns), sum(lns)/len(lns))
	print (total)
	return dict(contig_genes)
def gene_stats(deNovo, CpGs):
	count_in = 0
	count_out = 0
	total = 0
	for ctg in CpGs.keys():
		for pos in CpGs[ctg]:
			within = False
			total += 1
			for coord in deNovo.get(ctg, []):
				if pos in range(coord[0], coord[1] + 1):
					count_in += 1
					within = True
					break
			if not within:
				count_out += 1
	print (count_in, count_out, total)

def mapMeth(freq):
	meth_pos = defaultdict(list)
	with open(freq, "r") as f1:
		print (next(f1))
		for line in f1:
			data = line.split()
			ctg = data[0]
			start = int(data[1])
			end = int(data[2])
			freq = float(data[6])
			#if freq >= 0.7:
			#	meth_pos[ctg].append([start, end])
			meth_pos[ctg].append([start, end, freq]) 
	return(meth_pos)
def GC_percent(seq):
	# Initialize a count for GC base pairs
	gc_count = 0
	pattern = r"(CG)"
	matches = [match.start() for match in re.finditer(pattern, seq)]
	cg_num = len(matches)
	count_c = 0
	count_g = 0
	# Iterate through the DNA sequence
	for base in seq:
		if base == "G": count_g += 1
		if base == "C": count_c += 1
		if base == "G" or base == "C":
			gc_count += 1
	# Print the GC count
	percentage = (gc_count/len(seq))*100
	#print ("GC COUNT", count_c, count_g)
	cg2gc = None
	if count_c == 0:
		print ("No C!")
	if count_g == 0:
		print ("No G!")
	try:
		cg2gc = (cg_num/len(seq))/((count_c/len(seq))*(count_g/len(seq)))
	except:
		pass
	return(("{:.2f}%".format(percentage), "{:.2f}".format(cg2gc)))
			
def GC(gmark, assembly):
	ctg_seq = {}
	in_genes = ""
	out_genes = ""
	all_genome = ""
	for rec in SeqIO.parse(assembly, "fasta"):
		ctg_seq[rec.id] = str(rec.seq).upper()
		all_genome += str(rec.seq).upper()
	init_out = 0
	for ctg in gmark.keys():
		if ctg_seq.get(ctg, None):
			for gene in sorted(gmark[ctg], key=lambda tup: tup[0]):
				start = gene[0] - 1
				end = gene[1] + 1
				in_genes += ctg_seq[ctg][start:end]
				out_genes += ctg_seq[ctg][init_out:start - 1]
				init_out = end
	#print ("genes:", len(in_genes), "out_genes:", len(out_genes), "total", len(in_genes) + len(out_genes))
	cds = GC_percent(in_genes)
	non_cds = GC_percent(out_genes)
	genome = GC_percent(all_genome)
	return {"genome":genome, "CDS":cds, "Non-CDS":non_cds}

def metaeuk_anno(unipreds, gmark):
	total_gmark = 0
	for ctg in gmark.keys():
		total_gmark += len(gmark[ctg])
	total_meta = 0
	matching = 0
	for rec in SeqIO.parse(unipreds, "fasta"):
		data = rec.id.split("|")
		total_meta += 1
		ctg = data[1]
		uni_acc = data[0]
		start = int(data[6])
		end = int(data[7])
		for coord in gmark.get(ctg, []):
			if start in range(coord[0], coord[1] + 1) or end in range(coord[0], coord[1] + 1):
				matching += 1
				break
	print ("genemark: ", total_gmark, "metaeuk: ", total_meta, "matching: ", (matching/total_meta)*100)

def meth_annotation(mapped_meth, metaeuk_codon):
	#db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel/')
	cnt_go = 0
	cnt = 0
	anno_types = defaultdict(lambda: defaultdict(int))
	meth_genes = []
	no_meth_genes = []
	for rec in SeqIO.parse(metaeuk_codon, "fasta"):
		if len(rec.seq) >= 500:
			meth_flag = False
			data = rec.id.split("|")
			ctg = data[1]
			uni_acc = data[0]
			start = int(data[6])
			end = int(data[7])
			for coord in mapped_meth.get(ctg, []):
				if coord[0] in range(start, end + 1) and coord[1] in range(start, end + 1):
					meth_flag = True
			o_e = GC_percent(str(rec.seq).upper(), meth_flag)
			if o_e:
				if meth_flag:
					meth_genes.append(o_e)
				else:
					no_meth_genes.append(o_e)
	print (len(meth_genes), len(no_meth_genes))
	print ("METH:", sum(meth_genes) / len(meth_genes), min(meth_genes), max(meth_genes))
	print ("NO METH:", sum(no_meth_genes) / len(no_meth_genes), min(no_meth_genes), max(no_meth_genes))
	#SeqIO.write(meth_genes, open("methylated_cds.fasta"))
	"""
		go = db.get(uni_acc.encode('utf-8'))
		if go is not None:
			res = go.decode('utf-8').split("|")
			for r in res:
				desc = r.split("#")[1]
				desc_data = desc.split(":")
				anno_types[desc_data[0]][desc_data[1]] += 1
			cnt_go += 1
		cnt += 1
	db.close()
	print (cnt, cnt_go, cnt_go/cnt*100)
	graph_lists = {}
	full_names = {"C":"Cellular Component", "F":"Molecular Function", "P":"Biological Process"}
	for t in anno_types.keys():
		graph_lists[full_names[t]] = list(anno_types[t].items())
		with open("/disk3/ALL_ASSEMBLIES/57/{}_go.txt".format(t), "w") as f_in:
			for desc in anno_types[t]:
				#f_in.write("{}\t{}\n".format(anno_types[t][desc], desc))
				pass
	#pickle.dump(graph_lists, open("/disk3/ALL_ASSEMBLIES/graph_lists_115.p", "wb"))
	"""

if __name__ == '__main__':
	#methPos = mapMeth("/disk3/ALL_ASSEMBLIES/58/flye_3000/methylation_frequency.tsv")
	#meth_annotation(methPos, "/disk3/ALL_ASSEMBLIES/58/flye_3000/uniprotPredsMeta.codon.fas")

	#gmark = genemark("/disk4/GENE_FINDING/sample_70/genemark.gtf")
	#GC(gmark, "/disk4/GENE_FINDING/sample_70/sample_70_masked.fa")
	#unipreds = "/disk3/ALL_ASSEMBLIES/70/uniprotPreds.fas"
	#metaeuk_anno(unipreds, gmark)
	#genome("/disk4/GENE_FINDING/sample_70/sample_70_masked.fa")
	#genome("/disk3/ALL_ASSEMBLIES/ORIGINAL_ASSEMBLY/Enteromyxium_leei.fna")

	print ("####################")
	#multiGenome()
	#CpGs = find_CpG("/disk4/GENE_FINDING/sample_70/sample_70_masked.fa")

	#gene_stats(gmark, CpGs)
	#gene_CpG("/disk3/ALL_ASSEMBLIES/"/disk3/ALL_ASSEMBLIES/nematostella_cds.fna"")




	gene_CpG_test("/disk3/ALL_ASSEMBLIES/nematostella_cds.fna")
	#gene_CpG_all()
	#all_sample_genes()
	
