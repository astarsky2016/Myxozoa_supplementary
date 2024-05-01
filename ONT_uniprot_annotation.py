from Bio import SeqIO
import sys, plyvel, pickle, scispacy, spacy
from scispacy.linking import EntityLinker
from collections import defaultdict
from transformers import AutoTokenizer, AutoModelForTokenClassification
def map_text():
	cnt = 0
	desc_cnt = 0
	db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel_desc/', create_if_missing=True)
	with open("/disk5/UNIPROT/uniprot_trembl.dat", "r") as dat_in:
		records = []
		for line in dat_in:
			if line.startswith("ID"):
				record = {"acc":None, "desc":[]}
				cc_break = False
			if line.startswith("AC"):
				data = line.split()
				ids_plural = [el.strip() for el in data[1].split(";") if el.strip()]
				record["acc"] = ids_plural[0]
			if line.startswith("DE"):
				data = line.split("   ")
				data_clean = data[1].split("{")[0].strip()
				record["desc"].append(data_clean)
			if line.startswith("CC   ----------"):
				cc_break = True
				continue
			if line.startswith("CC"):
				if not cc_break:
					data_clean = line[9:].strip()
					record["desc"].append(data_clean)
			if line.startswith("KW"):
				data = line.split("   ")
				record["desc"].append(data[1].split("{")[0].strip())
			if line.startswith("//"):
				if record:
					cnt += 1
					if record["desc"]:
						records.append((record["acc"], " ".join(record["desc"])))
						desc_cnt += 1
					if len(records) == 100000:
						with db.write_batch() as wb:
							for rec in records:
								wb.put(rec[0].encode(), rec[1].encode())
						print ("batch written!")
						records = []
				record = {}
		if records:
			with db.write_batch() as wb:
				for rec in records:
					wb.put(rec[0].encode(), rec[1].encode())
			print ("final batch written!")
	print ("total", cnt, "desc", desc_cnt)
	db.close()

def map_GO():
	cnt = 0
	go_cnt = 0
	db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel/', create_if_missing=True)
	with open("/disk5/UNIPROT/uniprot_trembl.dat", "r") as dat_in:
		records = []
		for line in dat_in:
			if line.startswith("ID"):
				record = {"acc":None, "go":[]}
			if line.startswith("AC"):
				data = line.split()
				ids_plural = [el.strip() for el in data[1].split(";") if el.strip()]
				record["acc"] = ids_plural[0]
				if len(ids_plural) > 1:
					print (ids_plural)
					print (line)
					sys.exit(1)
			if line.startswith("DR   GO"):
				data = [el.strip() for el in line.split(";")]
				go_id = data[1]
				desc = data[2]
				record["go"].append({"go_id":go_id, "desc":desc})
			if line.startswith("//"):
				if record:
					cnt += 1
					if record["go"]:
						records.append(record)
						go_cnt += 1
					if len(records) == 100000:
						with db.write_batch() as wb:
							for rec in records:
								go_string = ""
								for go in rec["go"]:
									go_string += "|{}#{}".format(go["go_id"], go["desc"])
								wb.put(rec["acc"].encode(), go_string[1:].encode())
						print ("batch written!")
						records = []
				record = {}
		if records:
			with db.write_batch() as wb:
				for rec in records:
					go_string = ""
					for go in rec["go"]:
						go_string += "|{}#{}".format(go["go_id"], go["desc"])
					wb.put(rec["acc"].encode(), go_string[1:].encode())
			print ("final batch written!")
	print ("total", cnt, "go", go_cnt)
	db.close()
def mapMeth():
	promoters = defaultdict(list)
	methylated_true = defaultdict(list)
	methylated_false = defaultdict(list)
	max_ctg = (None, 0)
	with open("/disk3/ALL_ASSEMBLIES/57/flye_3000/promoters.txt", "r") as f_prom:
		for line in f_prom:
			if line.startswith("contig"):
				ctg = line.split(",")[0]
				ctg_size = int(line.split()[1])
				if ctg_size > max_ctg[1]:
					max_ctg = (ctg, ctg_size)
			if "Highly likely prediction" in line:
				data = [el.strip() for el in line.split()]
				pos = int(data[0])
				promoters[ctg].append(pos)
	met = 0
	unmet = 0
	tot = 0
	with open("/disk3/ALL_ASSEMBLIES/115/methylation_frequency.tsv", "r") as f1:
		print (next(f1))
		for line in f1:
			data = line.split()
			ctg = data[0]
			start = int(data[1])
			end = int(data[2])
			pos = int((start + end)/2)
			freq = data[6]
			tot += 1
			if float(freq) <= 0.3:
				methylated_false[ctg].append((start, end))
				unmet += 1
			if float(freq) >= 0.7:
				methylated_true[ctg].append((start, end))
				met += 1
	print (tot, met/tot, unmet/tot)
	c = max_ctg[0]
	for el in [methylated_true, methylated_false]:
		cnt = 0
		cnt_tot = 0
		for p in promoters[c]:
			cnt_tot += 1
			for tup in el[c]:
				if tup[0] in range(p - 100, p + 100) and tup[1] in range(p - 100, p + 100):
					cnt += 1
		#print (cnt, cnt_tot, cnt/cnt_tot)
	return (methylated_true, methylated_false, promoters)
def NER_analyze():
	f_txt = open("/disk3/ALL_ASSEMBLIES/sample_58.txt", "w")
	#nlp = spacy.load("en_core_sci_sm")
	#nlp.add_pipe("scispacy_linker", config={"resolve_abbreviations": True, "linker_name": "go"})
	db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel_desc/')
	text = ""
	for rec in SeqIO.parse("/disk3/ALL_ASSEMBLIES/58/flye_3000/uniprotPredsMeta.fas", "fasta"):
		data = rec.id.split("|")
		ctg = data[1]
		uni_acc = data[0]
		desc = db.get(uni_acc.encode('utf-8'))
		if desc is not None:
			text += desc.decode('utf-8')
			text += " "
	print (len(text))
	f_txt.write(text)
	f_txt.close()
	"""
	#text = "The expression of the TP53 gene was analyzed, and the protein p53 was detected."
	model_name = "dmis-lab/biobert-v1.1"
	tokenizer = AutoTokenizer.from_pretrained(model_name)
	model = AutoModelForTokenClassification.from_pretrained(model_name)
	# Tokenize the text
	tokens = tokenizer.tokenize(text)
	# Get the model's predictions
	inputs = tokenizer(text, return_tensors="pt")
	outputs = model(**inputs)
	# Extract predicted labels
	predicted_labels = outputs.logits.argmax(dim=2)
	# Map the labels back to entity names
	entity_labels = [model.config.id2label[label.item()] for label in predicted_labels[0]]
	# Print the text with entity labels
	for token, label in zip(tokens, entity_labels):
		print(f"Token: {token}, Label: {label}")
	doc = nlp(text)
	all_ents = defaultdict(int)
	cnt = 0
	for ent in doc.ents:
		all_ents[ent.text] += 1
	sorted_ny_len = sorted(all_ents.items(), key=lambda tup: tup[1], reverse=True)
	print (sorted_ny_len[0:50])
	"""
def metaeuk_anno_meth(m_true, m_false, promoters):
	db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel/')
	cnt_go = 0
	cnt = 0
	met_genes = set()
	unmet_genes = set()
	anno_types_high = defaultdict(lambda: defaultdict(int))
	anno_types_low = defaultdict(lambda: defaultdict(int))
	methylation_pattern = defaultdict(list)
	for rec in SeqIO.parse("/disk3/ALL_ASSEMBLIES/115/predsResults.fas", "fasta"):
		met_high = False
		met_low = False
		data = rec.id.split("|")
		ctg = data[1]
		uni_acc = data[0]
		start = int(data[6])
		end = int(data[7])
		for coord in m_true.get(ctg, []):
			if coord[0] in range(start, end + 1) and coord[1] in range(start, end + 1):
				methylation_pattern[uni_acc].append(1)
				met_high = True
		for coord in m_false.get(ctg, []):
			if coord[0] in range(start, end + 1) and coord[1] in range(start, end + 1):
				methylation_pattern[uni_acc].append(0)
				met_low = True
		go = db.get(uni_acc.encode('utf-8'))
		if go is not None:
			res = go.decode('utf-8').split("|")
			for r in res:
				desc = r.split("#")[1]
				desc_data = desc.split(":")
				if met_high:
					anno_types_high[desc_data[0]][desc_data[1]] += 1
				if met_false and not met_high:
					anno_types_low[desc_data[0]][desc_data[1]] += 1
			cnt_go += 1
		cnt += 1
	db.close()
	print (cnt, cnt_go, cnt_go/cnt*100)
	graph_lists_high = {}
	graph_lists_low = {}
	full_names = {"C":"Cellular Component", "F":"Molecular Function", "P":"Biological Process"}
	for t in anno_types_high.keys():
		graph_lists_high[full_names[t]] = list(anno_types_high[t].items())
	pickle.dump(graph_lists_high, open("/disk3/ALL_ASSEMBLIES/graph_lists_115_high.p", "wb"))
	for t in anno_types_low.keys():
		graph_lists_low[full_names[t]] = list(anno_types_low[t].items())
	pickle.dump(graph_lists_low, open("/disk3/ALL_ASSEMBLIES/graph_lists_115_low.p", "wb"))
	print ("low:", len(graph_lists_low), "high:", len(graph_lists_high))
	met_genes = 0
	unmet_genes = 0
	for acc in methylation_pattern.keys():
		if sum(methylation_pattern[acc]) > 0:
			met_genes += 1
			print (methylation_pattern[acc])
		else:
			unmet_genes += 1
	print ("Gene counts:", met_genes, unmet_genes)
def metaeuk_anno(m_true, m_false, promoters):
	db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel/')
	cnt_go = 0
	cnt = 0
	anno_types = defaultdict(lambda: defaultdict(int))
	for rec in SeqIO.parse("/disk3/ALL_ASSEMBLIES/57/flye_3000/uniprotPredsMeta.fas", "fasta"):
		data = rec.id.split("|")
		ctg = data[1]
		uni_acc = data[0]
		start = int(data[6])
		end = int(data[7])
		for coord in m_true.get(ctg, []):
			if coord[0] in range(start, end + 1) and coord[1] in range(start, end + 1):
				print ("Tuu je MET")
		for coord in m_false.get(ctg, []):
			if coord[0] in range(start, end + 1) and coord[1] in range(start, end + 1):
				print ("Tuu je NONE")
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

if __name__ == '__main__':
	#map_text()
	#NER_analyze()
	#map_GO()
	met_true, met_false, promo = mapMeth()
	#metaeuk_anno(met_true, met_false, promo)
	metaeuk_anno_meth(met_true, met_false, promo)
