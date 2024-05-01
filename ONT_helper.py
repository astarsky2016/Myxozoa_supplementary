from Bio import SeqIO, Entrez
from collections import defaultdict
import pickle, sys, csv, subprocess, os, gzip
from ete3 import NCBITaxa
ncbi = NCBITaxa()
Entrez.email = "your_email@example.com"
def getReadsAbove():
	threshold = 10000
	fastq = "/disk5/RUN1/sample_58/myxo_mapped.fastq.gz"
	reads = "/disk5/RUN1/sample_58/myxo_mapped_10000.fastq"
	batch = []
	with gzip.open(fastq, "rt") as handle:
		with open(reads, "a") as f_out:
			for rec in SeqIO.parse(handle, "fastq"):
				if len(rec.seq) >= threshold:
					batch.append(rec)
				if len(batch) >= 100000:
					SeqIO.write(batch, f_out, "fastq")
					batch = []
			if batch:
				SeqIO.write(batch, f_out, "fastq")
				batch = []
	print ("filtered")
def calculate_composite_score(bitscore, percent_identity, alignment_length, max_alignment_length):
	weights = {'bitscore': 0.7,'percent_identity': 4,'alignment_length': 2}
	# Normalize the parameters
	normalized_bits = bitscore
	normalized_percent_identity = percent_identity #already in percent
	normalized_alignment_length = alignment_length / max_alignment_length  # Normalize by dividing with the maximum alignment length
	# Calculate the composite score
	composite_score = (weights['bitscore'] * normalized_bits) + (weights['percent_identity'] * normalized_percent_identity) + (weights['alignment_length'] * normalized_alignment_length)
	return composite_score

def reduceHits(ctg_hits, fastas, f_out):
	#summarizing animal composition based on hit position and e_val
	animal_composition = {}
	for el in ctg_hits["contig_9155"]:
		print (el["qstart"], el["qend"])
	sys.exit(1)
	for ctg in ctg_hits.keys():
		#sorting hits by start position in contig
		hits = sorted(ctg_hits[ctg], key=lambda x: x["qstart"])
		pos_specific = {}
		num = 0
		for hit in hits:
			if hit["qstart"] > hit["qend"]:
				start = hit["qlen"] + 1 - hit["qstart"]
				end = hit["qlen"] + 1 - hit["qend"]
				hit["qstart"] = start
				hit["qend"] = end
				print ("changed to: {} - {}".format(start, end))
		prior = [set(range(hits[0]["qstart"], hits[0]["qend"])), calculate_composite_score(hits[0]["bitscore"], 
		hits[0]["pident"], hits[0]["length"], hits[0]["slen"])]
		pos_specific[0] = hits[0]
		for hit in hits[1:]:
			med = int((hit["qstart"] + hit["qend"])/2) 
			score = calculate_composite_score(hit["bitscore"],hit["pident"], hit["length"], hit["slen"])
			if med in prior[0]:
				if score > prior[1]:
					pos_specific[num] = hit
					prior[1] = score
			else:
				if hit["qstart"] > max(prior[0]): #razmisli o micanju
					num += 1
					pos_specific[num] = hit
					prior = [set(range(hit["qstart"], hit["qend"])), calculate_composite_score(hit["bitscore"],hit["pident"],
						hit["length"], hit["slen"])]
		animal_composition[ctg] = pos_specific
	max_diversity = (None, 0)
	for ctg in animal_composition.keys():
		animals = ""
		diversity = set()
		for pos in animal_composition[ctg].keys():
			hit = animal_composition[ctg][pos]
			try:
				l = ncbi.get_lineage(hit["staxid"])
				names = ncbi.get_taxid_translator(l)
				ranks = ncbi.get_rank(l)
				inv_rank = dict([(e[1], e[0]) for e in ranks.items()])
				animals += "##{}".format(names[inv_rank["class"]])
				diversity.add(names[inv_rank["class"]])
			except:
				pass
				#print ("missing data")
		if "Actinopteri" in diversity and "Myxozoa" in diversity: 
			print (ctg, animals)
			if len(diversity) > max_diversity[1]:
				max_diversity = (ctg, len(diversity))
			print ("###################")
	print (max_diversity, len(fastas[max_diversity[0]].seq))
	with(open(f_out, "w")) as fast_out:
		SeqIO.write([fastas[max_diversity[0]]], f_out, "fasta")
	print (f_out, "written")

def maxContigs(assembly):
	max_ctg = (None, 0)
	head, tail = os.path.split(assembly)
	for rec in SeqIO.parse(assembly, "fasta"):
		if len(rec.seq) > max_ctg[1]:
			max_ctg = (rec, len(rec.seq))
	print (max_ctg[1])
	max_file = os.path.join(head, "max_contig.fasta")
	blast_res = os.path.join(head, "max_contig.tab")
	with open(max_file, "w") as f_out:
		SeqIO.write([max_ctg[0]], f_out, "fasta")
	blast_command = [
	'blastn',
	'-query', max_file,
	'-db', '/disk4/GENBANK/NT_euk/nt_euk',
	'-out', blast_res,
	'-evalue', '1e-25',
	'-max_target_seqs', '10',
	'-max_hsps', '1',
	'-outfmt', '6 qseqid sseqid pident length slen qstart qend sstart send evalue bitscore staxids qstrand'
	]
	subprocess.run(blast_command, check=True)
	hits = defaultdict(list)
	with open(blast_res, "r") as blst:
		for line in blst:
			data = [el.strip() for el in line.split()]
			tx = data[11].split(";")
			hits[data[0]].append({"pident":float(data[2]), "length":int(data[3]), "slen":int(data[4]), 
			"qstart":int(data[5]), "qend":int(data[6]), "bitscore":int(data[10]), "staxid":int(tx[0])})
	reduceHits(hits)

def startFromBlast(blst_res, assembly):
	fastas = {}
	head, tail = os.path.split(assembly)
	max_file = os.path.join(head, "max_contig.fasta")
	for rec in SeqIO.parse(assembly, "fasta"):
		fastas[rec.id] = rec
	hits = defaultdict(list)
	with open(blst_res, "r") as blst:
		for line in blst:
			data = [el.strip() for el in line.split()]
			tx = data[11].split(";")
			hits[data[0]].append({"pident":float(data[2]), "length":int(data[3]), "slen":int(data[4]), 
			"qstart":int(data[5]), "qend":int(data[6]), "bitscore":int(data[10]), "staxid":int(tx[0])})
	reduceHits(hits, fastas, max_file)

def startFromDiamond(dmnd_res, assembly):
	#~/diamond blastx --query polished_assembly.fa.masked --db /disk5/NR/nrDMND --evalue 1e-20 --max-target-seqs 1 /
	#--sensitive --outfmt 6 qseqid sseqid pident length slen qstart qend sstart send evalue bitscore staxids > chimeric_diamond.out
	fastas = {}
	head, tail = os.path.split(assembly)
	max_file = os.path.join(head, "max_contig.fasta")
	for rec in SeqIO.parse(assembly, "fasta"):
		fastas[rec.id] = rec
	hits = defaultdict(list)
	with open(dmnd_res, "r") as blst:
		for line in blst:
			data = [el.strip() for el in line.split()]
			try:
				tx = data[11].split(";")
			except:
				print (data)
				continue
			hits[data[0]].append({"pident":float(data[2]), "length":int(data[3]), "slen":int(data[4]), 
			"qstart":int(data[5]), "qend":int(data[6]), "bitscore":float(data[10]), "staxid":int(tx[0]), "qlen":len(fastas[data[0]])})
	reduceHits(hits, fastas, max_file)
def fastq_screen(fasta):
	head, tail = os.path.split(fasta)
	top100 = dict(zip(range(100), ["None"]*100))
	for rec in SeqIO.parse(fasta, "fasta"):
		if len(rec.seq) > min(list(top100.keys())):
			top100_sort = sorted(top100.items(), key=lambda tup: tup[0])
			del top100[top100_sort[0][0]]
			top100[len(rec.seq)] = rec
	with open(os.path.join(head, "long_reads.fasta"), "w") as f_out:
		SeqIO.write(list(top100.values()), f_out, "fasta")
	print (list(top100.keys()))
def getProteins(dmnd_res):
	head, tail = os.path.split(dmnd_res)
	ids = set()
	with open(dmnd_res, "r") as blst:
		for line in blst:
			data = [el.strip() for el in line.split()]
			try:
				tx = data[11].split(";")
			except:
				continue
			for t in tx:
				try:
					l = ncbi.get_lineage(int(t))
					if 186623 not in l:
						ids.add(data[1])
				except:
					print ("tax issue")
	print (len(ids))
	# Fetch the protein sequences from NCBI using the accession numbers
	handle = Entrez.efetch(db="protein", id=list(ids), rettype="fasta", retmode="text")
	records = SeqIO.parse(handle, "fasta")
	# Save the protein sequences to a file
	SeqIO.write(records, os.path.join(head, "protein_sequences.fasta"), "fasta")
	# Close the handle
	handle.close()


if __name__ == '__main__':
	getReadsAbove()
	#maxContigs("/disk6/TMP_PARA/sample_108_new/flye_meta/BLOBS2/filtered_assembly/polished_assembly.fa")
	#startFromBlast("/disk6/TMP_PARA/sample_108_new/flye_meta/BLOBS2/filtered_assembly/chimeric_diamond.out", "/disk6/TMP_PARA/sample_108_new/flye_meta/BLOBS2/filtered_assembly/polished_assembly.fa.masked")
	#startFromDiamond("/disk6/TMP_PARA/sample_108_new/flye_meta/BLOBS2/filtered_assembly/chimeric_diamond.out", "/disk6/TMP_PARA/sample_108_new/flye_meta/BLOBS2/filtered_assembly/polished_assembly.fa.masked")
	#fastq_screen("/disk2/PAUL_PARA/X204SC22062342-Z01-F001_05/raw_data/Ceratomyx108/20220715_1648_6B_PAK14653_9e5215b8/fastq_pass/all_pass.fasta")
	#getProteins("/disk6/TMP_PARA/sample_108_new/flye_meta/BLOBS2/filtered_assembly/chimeric_diamond.out")
