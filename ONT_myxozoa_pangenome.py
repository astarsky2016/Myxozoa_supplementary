from Bio import SeqIO
import sys, plyvel, pickle, scispacy, spacy, glob, os
from scispacy.linking import EntityLinker
from collections import defaultdict
from transformers import AutoTokenizer, AutoModelForTokenClassification
def NER_analyze():
	f_txt = open("/disk3/ALL_ASSEMBLIES/sample_70.txt", "w")
	#nlp = spacy.load("en_core_sci_sm")
	#nlp.add_pipe("scispacy_linker", config={"resolve_abbreviations": True, "linker_name": "go"})
	db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel_desc/')
	text = ""
	for rec in SeqIO.parse("/disk3/ALL_ASSEMBLIES/70/uniprotPreds.fas", "fasta"):
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

def metaeuk_anno():
	db = plyvel.DB('/disk5/UNIPROT/uniprot_plyvel/')
	sample_results = {}
	for metaeuk_anno in glob.glob("/disk3/ALL_ASSEMBLIES/ORIGINAL_ASSEMBLY/*metaeuk.fas"):
		base = os.path.basename(metaeuk_anno)
		print (base)
		sample_name = base.replace("_metaeuk.fas", "")
		anno_types = defaultdict(lambda: defaultdict(int))
		for rec in SeqIO.parse(metaeuk_anno, "fasta"):
			data = rec.id.split("|")
			ctg = data[1]
			uni_acc = data[0]
			go = db.get(uni_acc.encode('utf-8'))
			if go is not None:
				res = go.decode('utf-8').split("|")
				for r in res:
					desc = r.split("#")[1]
					desc_data = desc.split(":")
					anno_types[desc_data[0]][desc_data[1]] += 1
		sample_results[sample_name] = anno_types
	db.close()
	full_names = {"C":"Cellular Component", "F":"Molecular Function", "P":"Biological Process"}
	function_dict = defaultdict(dict)
	counts_dict = defaultdict(dict)
	for samp in sample_results.keys():
		for category in sample_results[samp].keys():
			functions = set(sample_results[samp][category].keys())
			counts = sorted(list(sample_results[samp][category].items()), key=lambda tup: tup[1], reverse=True)
			function_dict[samp][category] = functions
			counts_dict[samp][category] = counts
		print (samp)
	#pickle.dump(dict(function_dict), open("/disk3/ALL_ASSEMBLIES/venn_data.p", "wb"))
	pickle.dump(dict(counts_dict), open("/disk3/ALL_ASSEMBLIES/counts_data.p", "wb"))
if __name__ == '__main__':
	NER_analyze()
	#metaeuk_anno()