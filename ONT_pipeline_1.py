from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
import os, glob, random, itertools, sys, subprocess, csv, pickle
from subprocess import PIPE, run
from collections import defaultdict
from Bio import SeqIO
import pandas as pd
def meth_format():
	f_out = open("/disk3/ALL_ASSEMBLIES/58/flye_3000/methylation_calls_57_tmp.tsv", "w")
	with open("/disk3/ALL_ASSEMBLIES/58/flye_3000/methylation_calls_57.tsv", "r") as f_in:
		for line in f_in:
			data = line.split()
			if len(data) == 11:
				f_out.write(line)
			else:
				break
	f_out.close()
def methPerContigs(freq, ctgs):
	head, tail = os.path.split(ctgs)
	all_ctg_len = {}
	for rec in SeqIO.parse(ctgs, "fasta"):
		all_ctg_len[rec.id] = len(rec.seq)
	#with open(os.path.join(head, "selected.faa"), "w") as selected_in:
	#	SeqIO.write(all_ctg, selected_in, "fasta")
	cnts_met = defaultdict(int)
	cnts_no = defaultdict(int)
	with open(freq, "r") as f_in:
		print(next(f_in))
		for line in f_in:
			data = line.split()
			ctg = data[0]
			start = data[1]
			end = data[2]
			fr = float(data[6])
			if fr >= 0.7:
				cnts_met[ctg] += 1
			if fr <= 0.1:
				cnts_no[ctg] += 1
	tot_met = sum(cnts_met.values())
	tot_no = sum(cnts_no.values()) 
	avg = tot_met/(tot_met + tot_no)
	all_ctg = []
	for ctg in cnts_met.keys():
		avg_ctg = cnts_met[ctg]/(cnts_met[ctg] + cnts_no[ctg])
		all_ctg.append((abs(avg_ctg - avg), ctg, avg_ctg*100))
	sorted_ctg = sorted(all_ctg, key=lambda tup: tup[0], reverse=True)
	print (sorted_ctg[0:2])
	print (avg*100)
	print (tot_met, cnts_met[sorted_ctg[0][1]])
	sorted_contigs = sorted(cnts_met.items(), key=lambda tup: tup[1], reverse=True)
	print (sorted_contigs[0:5], [all_ctg_len[el[0]] for el in sorted_contigs[0:5]])
def methSanity():
	cnt = 0
	cnt_met = 0
	cnt_hypo = 0
	cnt_med = 0
	cov = []
	cov_med = []
	cov_hyp = []
	cov_tot = []
	hints = []
	#hints_108 = pickle.load(open("108_called.p", "rb"))
	with open("/disk6/TMP_PARA/sample_115/assembly5000/methylation_frequency.tsv", "r") as f_in:
		print(next(f_in))
		for line in f_in:
			data = line.split()
			cnt += 1
			if data[6] == "1.000":
				hints.append((data[0], data[1], data[2]))
			if data[5] == data[4]:
				cov_tot.append(int(data[4]))
			if float(data[6]) >= 0.8:
				cnt_met += 1
				cov.append(int(data[4]))
			else:
				if float(data[6]) <= 0.2:
					cnt_hypo += 1
					cov_hyp.append(int(data[4]))
				else:
					if float(data[6]) >= 0.3:
						cnt_med += 1
						cov_med.append(int(data[4]))
	print (cnt_met, cnt_med, cnt_hypo)
	print ((cnt_met/cnt)*100)
	print ((cnt_hypo/cnt)*100)
	print ("full", sum(cov)/len(cov), len(cov))
	print ("med", sum(cov_med)/len(cov_med), len(cov_med))
	print ("none", sum(cov_hyp)/len(cov_hyp), len(cov_hyp))
	print (len(cov_tot), sum(cov_tot)/len(cov_tot), max(cov_tot), max(cov))
	print ("###################")
	#print(len(hints), len(hints_108), len(set(hints).intersection(set(hints_108))))
	#pickle.dump(hints, open("108_called.p", "wb"))
def reformatIndex(index_file):
	new_idx = index_file.replace(".readdb", ".readdb_new")
	with open(new_idx, "w") as idx_in:
		with open(index_file, "r") as idx:
			for line in idx:
				new_line = line.replace("fast5_pass/", "/disk6/TMP_PARA/sample_57_data/fast5_pass")
				idx_in.write(new_line)
	print ("done!")
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
def filterOnly(fastas, workdir):
	threshold = 3000
	filtered = os.path.join(workdir, "stage1_over3000.fasta")
	batch = []
	with open(filtered, "a") as f_out:
		for rec in SeqIO.parse(fastas, "fasta"):
			if len(str(rec.seq)) >= threshold:
				batch.append(rec)
			if len(batch) >= 100000:
				SeqIO.write(batch, f_out, "fasta")
				batch = []
		if batch:
			SeqIO.write(batch, f_out, "fasta")
			batch = []
	print ("4. length filtered")
def firstStageAssembly(fastq_reads, workdir):
	threshold = 3000
	out_dir = os.path.join(workdir, "meta_assembly")
	os.mkdir(out_dir)
	filtered = os.path.join(workdir, "stage1_reads_over3000.fasta")
	batch = []
	with open(filtered, "a") as f_out:
		for rec in SeqIO.parse(fastq_reads, "fasta"):
			if len(rec.seq) >= threshold:
				batch.append(rec)
			if len(batch) >= 100000:
				SeqIO.write(batch, f_out, "fasta")
				batch = []
		if batch:
			SeqIO.write(batch, f_out, "fasta")
			batch = []
	print ("4. length filtered")
	assemby_cmd = ["flye", "--nano-hq", filtered, "--out-dir", out_dir, "--threads", "32"]
	assemble = run(assemby_cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("1. ASSEMBLED")
	return (os.path.join(out_dir, "assembly.fasta"), out_dir)
def fishFIltering(reads, workDir):
	# Run minimap2 command for fish mapping of contigs
	mapping = os.path.join(workDir, "fish_mapped_sorted.bam")
	tmpfile = os.path.join(workDir, "mapping.tmp")
	command = "minimap2 -ax map-ont -t 16 -I 50G /disk4/filter_fish/fish_filter.fna {} | samtools sort -o {} -T {}".format(reads, mapping, tmpfile)
	os.system(command)
	print ("1. mapped")
	# Run samtools index command
	samtools_index_command = ["samtools", "index", mapping, "-@", "12"]
	run(samtools_index_command, check=True)
	print ("2. indexed")
	#1. stage filtering - unmapped reads 
	fish_unmapped_reads_bam = os.path.join(workDir, "fish_unmapped_reads.bam")
	view_command = ["samtools","view","-b","-f", "4",mapping]
	with open(fish_unmapped_reads_bam, "wb") as output_file:
		subprocess.run(view_command, check=True, stdout=output_file)
	# Run samtools fastq command
	fish_unmapped_reads_fastq = fish_unmapped_reads_bam.replace(".bam", ".fasta")
	fasta_command = ["samtools","fasta",fish_unmapped_reads_bam]
	with open(fish_unmapped_reads_fastq, "w") as output_file:
		subprocess.run(fasta_command, check=True, stdout=output_file)
	print ("3. first stage fastq created")
	return fish_unmapped_reads_fastq
def consensus(task):
	tag = task[0]
	contigs = " ".join(task[1])
	bam = task[2]
	hdf = task[3]
	batch = 200
	threads = 2
	cmd = "medaka consensus {} {} --model r941_prom_high_g4011 --batch 200 --threads 2 --region {}".format(bam, hdf, contigs)
	os.system(cmd)
	return(len(task[1]))
def medaka_par(assembly, basecalls, workDir):
	hdf_dir = os.path.join(workDir, "HDF")
	contigs = []
	total = 0
	limit = 0
	for rec in SeqIO.parse(assembly, "fasta"):
		total += 1
		if len(rec.seq) >= limit:
			max_ctg = len(rec.seq)
			limit = len(rec.seq)
		contigs.append(rec)
	print (len(contigs), "max:", max_ctg)
	os.mkdir(hdf_dir)
	mapping = os.path.join(workDir, "calls_to_draft")
	print ("mapping")
	cmd = ["mini_align", "-i", basecalls, "-r", assembly, "-P", "-m", "-p", mapping, "-t", "36"]
	result = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("reads mapped")
	tasks = list(chunks([c.id for c in contigs], 10))
	chunk_elements = [(pos, el, mapping + ".bam", os.path.join(hdf_dir, "contigs_{}.hdf".format(pos))) for pos, el in enumerate(tasks)]
	print (len(chunk_elements))
	with mp.Pool(20) as pool:
		for result in pool.map(consensus, chunk_elements):
			print ("finished chunks", result)
	polished_assembly = os.path.join(workDir, "polished_assembly_clean.fa")
	stitch_cmd = "medaka stitch --threads 20 {}/*.hdf {} {}".format(hdf_dir, assembly, polished_assembly)
	os.system(stitch_cmd)
	print ("7. ASSEMBLY POLISHED")
	return polished_assembly
def blob_prep(assembly):
	# Run samtools view command to get fish unmapped contigs
	head, tail = os.path.split(assembly)
	"""
	fish_unmapped_contigs_bam = os.path.join(head, "fish_unmapped_contigs.bam")
	mapping = os.path.join(head, "mapped_reads_stage1.bam")
	view_command = ["samtools","view","-b","-f", "4",sorted_bam]
	with open(fish_unmapped_contigs_bam, "wb") as output_file:
		subprocess.run(view_command, check=True, stdout=output_file)
	# Run samtools fasta command
	fish_unmapped_contigs_fasta = fish_unmapped_contigs_bam.replace(".bam", ".fasta")
	fasta_command = ["samtools","fasta",fish_unmapped_contigs_bam]
	with open(fish_unmapped_contigs_fasta, "w") as output_file:
		subprocess.run(fasta_command, check=True, stdout=output_file)
	print ("fasta created")
	# re-align original reads to unmapped contigs for coverage
	command = "minimap2 -ax asm20 -t 16 -I 50G {} {} | samtools sort -o {} -T {}".format(filtered_stage1, original_reads, mapping, tmpfile)
	run(command, shell=True, check=True)
	print ("reads re-aligned")
	"""
	# diamond blastx command
	print ("diamond")
	diamond_out = os.path.join(head, "polished_diamond.tab")
	diamond_command = ["/home/astar/diamond","blastx","-q", assembly,"-d", "/disk6/TSA_new/tsa_diamond.dmnd",
	"--outfmt", "6 qseqid staxids bitscore qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
	"--sensitive","--max-target-seqs", "1","--evalue", "1e-25","--threads", "24","-o", diamond_out]
	cmd_str = " ".join(diamond_command)
	os.system(cmd_str)
	print ("diamond success")
	# blastn command
	print ("blast")
	blastn_out = os.path.join(head, "polished_blastn.tab")
	blastn_command = ["blastn","-db", "/disk4/GENBANK/NT_euk/nt_euk",
	"-query", assembly,"-outfmt", '"6 qseqid staxids bitscore std"',
	"-max_target_seqs", "10","-max_hsps", "1","-evalue", "1e-25",
	"-num_threads", "32","-out", blastn_out]
	cmd_str = " ".join(blastn_command)
	os.system(cmd_str)
	print ("blastn success")
	print ("3. hits assigned")
def assemble_final(reads):
	head, tail = os.path.split(reads)
	out_dir = os.path.join(head, "assembly_2")
	assembly_cmd = ["flye", "--nano-hq", reads, "--out-dir", out_dir, "--scaffold", "--threads", "24"]
	assemble = run(assembly_cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("5. ASSEMBLED SINGLE")
def makeBed(bed_file, contigs, mapping):
	"""
	with open(bed_file, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', lineterminator='\n')
		for ctg in contigs:
			writer.writerow([ctg[0], "0", "{}".format(ctg[1] - 1)])
	"""
	mapped_reads = bed_file.replace(".bed", "_mapped.bam")
	mapped_reads_fastq = mapped_reads.replace(".bam", ".fastq")
	cmd = "samtools view -b -F 4 -L {} {} -o {}".format(bed_file, mapping, mapped_reads)
	cmd2 = "samtools fastq {} > {}".format(mapped_reads, mapped_reads_fastq)
	try:
		#subprocess.run(cmd, shell=True, check=True)
		print ("reads extracted")
		#subprocess.run(cmd2, shell=True, check=True)
		print ("reads recorded")
		assemble_2 = mapped_reads_fastq.replace(".fastq", "2000.fastq")
		"""
		with open(assemble_2, "a") as out:
			batch2000 = []
			for rec in SeqIO.parse(mapped_reads_fastq, "fastq"):
				if len(rec.seq) >= 2000:
					batch2000.append(rec)
					if len(batch2000) >= 100000:
						SeqIO.write(batch2000, out, "fastq")
						batch2000 = []
			if batch2000:
				SeqIO.write(batch2000, out, "fastq")
		print ("assemby 2 dataset prepared")
		"""
		print(assemble_2)
		assemble_final(assemble_2)
	except subprocess.CalledProcessError as e:
		print(f"Error executing command: {e}")
def extractFilteredStage2(csvFile, mapping, assembly):
	# Read the CSV file into a DataFrame
	df = pd.read_csv(csvFile)
	# Print the DataFrame
	contigs = df["id"].tolist()
	nonzero_coverage = df.loc[df["reads-polished.sorted_cov"] != 0, ["id", "length"]]
	nonzero_coverage_contigs = nonzero_coverage.values.tolist()
	print (len(contigs), len(nonzero_coverage_contigs))
	print (nonzero_coverage_contigs[0:5])
	selection = [el[0] for el in nonzero_coverage_contigs]
	clean = []
	for rec in SeqIO.parse(assembly, "fasta"):
		if rec.id in selection:
			if len(rec.seq) >= 2000:
				clean.append(rec)
	with open(assembly.replace(".fa", "_stage2.fa"), "w") as out:
		SeqIO.write(clean, out, "fasta")
	#makeBed(csvFile.replace(".csv", ".bed"), nonzero_coverage_contigs, mapping)
def comparePolished(a1, a2):
	s1 = set()
	s2 = set()
	tot1 = 0
	for rec in SeqIO.parse(a1, "fasta"):
		s1.add(rec.id)
		tot1 += len(rec.seq)
	tot2 = 0
	for rec in SeqIO.parse(a2, "fasta"):
		s2.add(rec.id)
		tot2 += len(rec.seq)
	print (len(s1))
	print (len(s2))
	print (len(s1.intersection(s2)))
	print (tot1, tot2)
def BUSCO_compare():
	buscos = defaultdict(set)
	single_busc = defaultdict(int)
	for f in glob.glob("/disk3/ALL_ASSEMBLIES/ORIGINAL_ASSEMBLY/CORE_BUSCOS/*.tsv"):
		name = os.path.basename(f).split(".")[0]
		with open(f, "r") as f_in:
			for line in f_in:
				if line.startswith("#"):
					continue
				else:
					data = [el.strip() for el in line.split()]
					if "Complete" in data:
						buscos[name].add(data[0])
						single_busc[data[0]] += 1
					if "Duplicated" in data:
						buscos[name].add(data[0])
						single_busc[data[0]] += 1
	all_buscos = set.union(*buscos.values())
	print (len(set.intersection(*buscos.values())))
	print (len(set.union(*buscos.values())))
	print (len(buscos["sample_70"]))
	test = all_buscos.difference(buscos["sample_70"])
	for samp in buscos.keys():
		print (samp, len(buscos[samp].intersection(buscos["sample_70"]))/len(buscos[samp]))
	sys.exit(1)
	l1 = []
	l2 = []
	for n in buscos.keys():
		if "sample" in n:
			print (n)
			l1.append(buscos[n])
		else:
			l2.append(buscos[n])
	u1 = set.intersection(*l1)
	u2 = set.intersection(*l2)
	print ("intersection", len(u1), len(u2))
	u1 = set.union(*l1)
	u2 = set.union(*l2)
	print ("union", len(u1), len(u2))
	print (len(buscos["sample_108"]), len(buscos["sample_115"]), len(buscos["sample_108"].intersection(buscos["sample_115"])))
	sorted_single = sorted(single_busc.items(), key=lambda tup: tup[1], reverse=True)
	print (sorted_single[0:5])
	print(sorted_single[-5:])
	subs = []
	for el in buscos["sample_70"]:
		subs.append(single_busc[el])
	print (max(subs), min(subs))
	print (sorted(subs))
if __name__ == '__main__':
	#methPerContigs("/disk6/TMP_PARA/sample_70/flye9000/methylation_frequency.tsv", "/disk3/ALL_ASSEMBLIES/70/sample_70_ceratomyxa_sp.fna")
	#methSanity()
	#sys.exit(1)
	#reformatIndex("/disk6/TMP_PARA/sample_57_data/all_pass.fastq.index.readdb")
	#meth_format()
	#sys.exit(1)
	reads = "/disk5/RUN1/sample_58/all_pass.fastq"
	workDir = "/disk5/RUN1/70_mapping/578"
	#fastas = "/disk6/TMP_PARA/sample_70/fish_unmapped_reads.fasta"
	#fastq_reads = fishFIltering(reads, workDir)
	#fastq_reads = "/disk5/RUN1/sample_58/NEW/fish_unmapped_reads.fastq"
	#res = firstStageAssembly(fastq_reads, workDir)
	#fastq_reads = "/disk5/RUN1/sample_57/sample_115_guided/ref_mapped_reads.fastq"
	#filterOnly(reads, workDir)
	#print (fastq_reads)
	contigs = "/disk3/ALL_ASSEMBLIES/58/flye_3000/assembly.fasta"
	out_dir = "/disk3/ALL_ASSEMBLIES/58/flye_3000"
	#out_dir = "/disk6/TMP_PARA/sample_70/canu_corr_all/stage1_trim"
	#contigs, out_dir = firstStageAssembly(fastq_reads, workDir)
	#polished = medaka_par(contigs, reads, out_dir)
	#for p in ["/disk5/RUN1/sample_57/NEW/meta_assembly/58_115mapped_contigs.fasta", "/disk5/RUN1/sample_58/NEW/meta_assembly/115mapped_contigs.fasta"]:
	#	blob_prep(p)
	a1 = "/disk6/TMP_PARA/sample_115/assembly5000/polished_assembly_clean.fa"
	a2 = "/disk6/TMP_PARA/sample_115/assembly5000/assemblystage2_blob_filt.fasta"
	#comparePolished(a1, a2)
	actinopteri_filt = "/disk3/ALL_ASSEMBLIES/58/flye_3000/sample_58_actinopteri.csv"
	mapping = "/disk6/TMP_PARA/sample_115/assembly5000/mapped_reads_stage1.bam"
	assembly = "/disk3/ALL_ASSEMBLIES/58/flye_3000/polished_assembly_clean.fa"
	#extractFilteredStage2(actinopteri_filt, mapping, assembly)
	BUSCO_compare()
