from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
import os, glob, random, itertools, sys, subprocess
from subprocess import PIPE, run
from collections import defaultdict
def filt():
	above = []
	for rec in SeqIO.parse("/disk6/TMP_PARA/sample_57_new/flye_meta/SANITY/actinopteri_filter_mapped.fastq", "fastq"):
		if len(rec.seq) >= 3000:
			above.append(rec)
	with open("/disk6/TMP_PARA/sample_57_new/flye_meta/SANITY/actinopteri_filter_mapped_3000.fasta", "w") as f_out:
		SeqIO.write(above, f_out, "fasta")
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]
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
if __name__ == '__main__':
	#filt()
	basecalls = "/disk2/PAUL_PARA/X204SC22062342-Z01-F001_01/raw_data/Ellipsomyx70/20220715_1647_6A_PAK11434_5c0f67ef/all_pass.fastq"
	workDir = "/disk6/TMP_PARA/sample_70/canu_corr_all/stage1_trim"
	medaka_par("/disk6/TMP_PARA/sample_70/canu_corr_all/stage1_trim/assembly.fasta", basecalls, workDir)
