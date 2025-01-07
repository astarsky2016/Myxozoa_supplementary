# Myxozoa_supplementary
Supplementary Python code for Myxozoa ONT data processing

Oxford-Nanopore MinION sequencing has been used to produce complete and near-complete assemblies of five myxozoan genomes free from host DNA contamination. It was also used to perform DNA methylation analysis indicated the presence of DNA methylation in the GC-rich regions within gene bodies. Here are main helper scripts that have been written and used for this tasks. This project contains main Python scripts used to analyze ONT fast5 datasets and produce results together with Jupyter notebook files used to generate all the figures in the manuscript. The methods and logic behind this work is described in the manuscript. These scripts are only working excerpts not meant to be used outside of the scope of work described in the manuscript. If one wishes to reproduce the work, please change input filenames manually, since these scripts contain hardcoded file names.

DNMT_annotation.ipynb - Jupyter notebook used for parsing of Panther and MetaEuk annotations of Myxozoa genome assemblies with the goal of ientifying BER, TET and DNMT parts of methylation machinery

GO_chart.ipynb - Jupyter notebook used to plot most important GO categories identified in Myxozoa genome assemblies

Gene_CpG_all.ipynb - CpG analysis of Myxozoa genome assemblies including KDE and bimodal distribution visualization

KDE_method.ipynb - Jupyter notebook that has been used to visualize observed bimodal distribution of CpG O/E values, includes KDE and Venn diagrams of GO distribution according to methylated/unmethylated parts of the genome

ONT_Silva_taxonomy.py - Helper script used in the process of analizing 18S rRNA content of samples based on SILVA ribosomal RNA database (locally downloaded)

ONT_helper.py - Helper script that has been used for ONT long read assessment prior to assembly. Long reads have been searched using Blastn and Diamond blast and the resulting output has been analyzed for the presence of genes that exhibit homology to known Myxozoa genes

ONT_methylation_distribution.py - GC content calculation, CpG O/E gene distribution across entire length of assembled Myxozoa genomes

ONT_myxozoa_pangenome.py - Named Entity Recognition based analysis of Myxozoa genome annotations

ONT_pipeline_1.py - Script that has been used in various stages of genome analyses: assembling of fish filtered reads, parallelization of Medaka polishing (additional round) of Fyle genome assemblies, mapping of methylation calling results obtained by nanopolish onto assembled genomes and preparing data for BlobTools additional round of filtering

ONT_pipeline_2.py - Medaka ploshing step parallelization script

ONT_prepare_data.py - Helper script used for alignemnt of ONT reads to fish genomes - later replaced by ONT_pipeline_1 script

ONT_uniprot_annotation.py - Script used to map Uniprot and GO based annotations onto assembled genomes with tagged methylation sites

all_sample_GO.ipynb - Jupyter notebook used to create radar plots and graphs summarizing GO categories per assembled genome

gene_CpG_test.ipynb - Jupyter notebook used in some analyses looking into distribution of CpG sites across annotated genes and looking itno gene lenght in particular

housekeeping.ipynb - NER analysis and Wordcloud depiction of some housekeeping genes

keywords.ipynb - KeyBERT based keyword analysis of PANTHER based genome annotations

methylation_signature.ipynb - Creation of heatmaps based on methylation frequancy data obtained using nanopolish

selected_GO.ipynb - Creation of interactive html files with scatterplots of identified GO

selected_GO_sideByside.ipynb - Creation of interactive html files with side-by-side scatterplots of identified GO for easier comparison

stacked_meth.ipynb - Creation of stacked barcharts displaying numbers of methylated/unmethylated genes per assembled Myxozoa genome

venn_myxo.ipynb - Creating Venn and SuperVenn diagrams showing shared GO across assembled Myxozoa genomes 
