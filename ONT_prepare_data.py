from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
import os, glob, random, itertools, sys, subprocess
from subprocess import PIPE, run
from collections import defaultdict

def alignToHost(basecalled_fastq, fast5_dir, reference_genomes, seq_summary):
	head, tail = os.path.split(basecalled_fastq)
	porechopped_fasta = os.path.join(head, "porechopped_reads.fasta")
	sorted_reads = os.path.join(head, "reads-ref.sorted.bam")
	unmapped_reads = os.path.join(head, "unmapped_reads.bam")
	unmapped_fasta = os.path.join(head, "unmapped_reads.fasta")
	tmp_file = os.path.join(head, "reads.tmp")
	"""
	nano_idx = ["/disk2/PAUL_PARA/nanopolish/nanopolish", "index", "-s", seq_summary, "-d", fast5_dir, basecalled_fastq]
	indexing_reads = run(nano_idx, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("1. READS INDEXED")
	
	align_cmd = ["minimap2", "-ax", "map-ont", "-t", "16", "-I", "50G", reference_genomes, basecalled_fastq, "|", "samtools",
	"sort", "-o", sorted_reads, "-T", tmp_file]

	alignment = run(align_cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("2. ALIGNED")

	index_cmd = ["samtools", "index", sorted_reads]
	indexing = run(index_cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("3. INDEXED SORTED")
	
	extract_cmd = "samtools view -b -f 4 {} > {}".format(sorted_reads, unmapped_reads)
	res = os.system(extract_cmd)
	#extracting = run(extract_cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	#Tu ima problema sa izvršavanjem, prekine prerano
	fasta_cmd = "samtools fasta {} > {}".format(unmapped_reads, unmapped_fasta)
	extracting = os.system(fasta_cmd)
	print ("4. UNMAPPED EXTRACTED")
	"""
	#tu možda treba sa seqkit razbiti na više manjih fasta ako je original prevelik
	split_reads = os.path.join(head, "SPLIT_CHUNKS")
	reads_split = breakFasta(unmapped_fasta, split_reads)
	completed = parallel_porechop(split_reads)
	#pore_clean_cmd = ["porechop", "-i", unmapped_fasta, "-o", porechopped_fasta, "--discard_middle"]
	#clean_reads = run(pore_clean_cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("5. CLEANED")
def assembleParasite(porechopped_reads):
	head, tail = os.path.split(porechopped_reads)
	out_dir = os.path.join(head, "assembly")
	os.mkdir(out_dir)
	assemby_cmd = ["flye", "--nano-hq", reads, "--out-dir", out_dir, "--threads", "24"]
	assemble = run(assemby_cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("6. ASSEMBLED")
	return out_dir
def metaEuk(polished_assembly):
	nr_DB = "/disk4/GENBANK/nr.fnaDB"
	head, tail = os.path.split(polished_assembly)
	polishedDB = os.path.join(head, "polishedDB")
	metaeuk_res = os.path.join(head, "nrPredsMeta")
	tmp_meta = os.path.join(head, "tmp_metaeuk")
	cmd1 = ["mmseqs", "createdb", polished_assembly, polishedDB, "--dbtype", "2"]
	assemble = run(cmd1, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("8. ASSEMBLED DB CREATED")
	cmd2 = ["metaeuk", "easy-predict", polishedDB, nr_DB, metaeuk_res, tmp_meta]
	print ("9. METAEUK PROTEINS ASSIGNED")
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
def medaka_par(assembly, basecalls, limit, workDir):
	hdf_dir = os.path.join(workDir, "HDF")
	contigs = []
	total = 0
	for rec in SeqIO.parse(assembly, "fasta"):
		total += 1
		if len(rec.seq) >= limit:
			contigs.append(rec)
	print (len(contigs), total)
	length_filtered = os.path.join(workDir, "contigs_over_{}.fasta".format(limit))
	os.mkdir(hdf_dir)
	mapping = os.path.join(workDir, "calls_to_draft")
	SeqIO.write(contigs, open(length_filtered, "w"), "fasta")
	print ("mapping")
	cmd = ["mini_align", "-i", basecalls, "-r", length_filtered, "-P", "-m", "-p", mapping, "-t", "48"]
	result = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("reads mapped")
	tasks = list(chunks([c.id for c in contigs], 10))
	chunk_elements = [(pos, el, mapping + ".bam", os.path.join(hdf_dir, "contigs_{}.hdf".format(pos))) for pos, el in enumerate(tasks)]
	print (len(chunk_elements))
	with mp.Pool(20) as pool:
		for result in pool.map(consensus, chunk_elements):
			print ("finished chunks", result)
	polished_assembly = os.path.join(workDir, "polished_assembly.fa")
	stitch_cmd = "medaka stitch --threads 20 {}/*.hdf {} {}".format(hdf_dir, assembly, polished_assembly)
	os.system(stitch_cmd)
	print ("7. ASSEMBLY POLISHED")
	return polished_assembly
def panther(fasta):
	tail, head = os.path.split(fasta)
	outfile = os.path.join(tail, "panther.out")
	cmd = ["/disk4/PANTHER_v17/pantherScore2.2.pl", "-l", "/disk4/PANTHER_v17/PANTHER17.0/", 
	"-D", "B", "-V", "-H", "/disk2/hmmer-3.3.2/bin/", "-i", fasta, "-o", outfile, "-n", "-s"]
	result = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print ("done!")
	return(outfile)
def polish(f):
	tail, head = os.path.split(f)
	porechopped_fasta = os.path.join(tail, head.replace(".fasta", "_chopped.fasta"))
	cmd = ["porechop", "-i", f, "-o", porechopped_fasta, "--discard_middle"]
	result = run(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
	print (porechopped_fasta)
	return porechopped_fasta
def parallel_porechop(d):
	outfiles = []
	fastas = glob.glob(os.path.join(d, "*.fasta"))
	print ("running")
	with mp.Pool(10) as pool:
		for result in pool.map(polish, fastas):
			outfiles.append(result)
	print ("done!")
	return 1
def breakFasta(f,d):
	os.mkdir(d)
	cnt = 0
	batch_counter = 0
	limit = 1000000
	batch = []
	for rec in SeqIO.parse(f, "fasta"):
		batch.append(rec)
		cnt += 1
		if cnt == limit:
			SeqIO.write(batch, open(os.path.join(d, "part_{}.fasta".format(batch_counter)), "w"), "fasta")
			batch = []
			batch_counter += 1
			cnt = 0
			print (batch_counter)
	if batch:
		SeqIO.write(batch, open(os.path.join(d, "part_{}.fasta".format(batch_counter)), "w"), "fasta")
	print ("done!")
	return 1
if __name__ == '__main__':
	#1
	chopped_reads = alignToHost("/disk6/TMP_PARA/sample_58/all_pass.fastq", "/disk6/TMP_PARA/X204SC22062342-Z01-F001_02/raw_data/Ceratomyxa58/20220715_1647_2F_PAK14463_e1513935/fast5_pass", 
		"/disk6/HOST_REF/scia_ref/sciaenidae.fna", "/disk6/TMP_PARA/X204SC22062342-Z01-F001_02/raw_data/Ceratomyxa58/20220715_1647_2F_PAK14463_e1513935/sequencing_summary_PAK14463_09b04042.txt")
	#breakFasta("/disk6/TMP_PARA/sample_57/unmapped_reads.fasta", "/disk6/TMP_PARA/sample_57/fastq_chunks")
	#parallel_porechop()
	#2
	#assembly_dir = assembleParasite(chopped_reads)
	#3
	#polished_assembly = medaka_par(os.path.join(assembly_dir, "assembly.fasta"), chopped_reads, 5000, assembly_dir)
	#4
	#metaEuk(polished_assembly)
	#panther("/disk6/TMP_PARA/sample_58/assembly/candidate_58_parasite.faa")
