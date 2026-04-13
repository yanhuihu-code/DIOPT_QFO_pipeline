# DIOPT_QFO_pipeline
This pipeline was designed to process the ortholog prediction result from multiple algorithms that were submitted to QFO (Quest for ortholog) benchmarking site. Then integrates the results with the option to supplement additional orthologous relationships from manually curated resources, including HGNC, RGD, and ZFIN.  Then calculates the voting scores and annotates the voting scores, whether it is the best score with forward and reverse searches, and the confidence level.  The user needs to download the prediction files to the local computer, make the vocabulary files regarding 1.) the species to be included, 2.) the prediction method, and corresponding file names. 
The files needed to be downloaded and run the pipeline are as follows
1.	Ortholog prediction results from QFO benchmarking site https://orthology.benchmarkservice.org/proxy/
Next are the examples of files used for the DIOPT assembly for AGR from the most recent QFO release (2022): https://orthology.benchmarkservice.org/proxy/projects/2022/
1.)	Compara: https://b2share.eudat.eu/records/mnbxy-pnx93/files/qfo_2022_02.default.orthoxml.xml.gz?download=1
2.)	Hieranoid: https://b2share.eudat.eu/records/rx86s-23p86/files/pairs_hieranoid-diamond.tsv?download=1
3.)	InParanoid: https://b2share.eudat.eu/records/ca8wf-v9w65/files/InParanoid_QFO22.pairs_1.gz?download=1
4.)	OMA: https://b2share.eudat.eu/records/w3r7b-da890/files/OMA.2.5.0-VPairs.txt.gz?download=1
5.)	OrthoFinder: https://b2share.eudat.eu/records/d8nrr-9fs80/files/OrthoFinder.txt?download=1
6.)	Panther: https://b2share.eudat.eu/records/c9wp1-7dh03/files/PANTHER18_all_1?download=1
7.)	Phylome: https://b2share.eudat.eu/records/yvdbk-8te46/files/orthologs.farthest.30Cl.2022.txt.gz?download=1
8.)	SonicParanoid: https://b2share.eudat.eu/records/e5k5r-vsk83/files/qfo22-challenge.sp2-sens.tsv.gz?download=1

2.	Additional ortholog predictions from manually curated resources (optional)
1.)	HGNC (Human-mouse): use custom download page at HGNC and choose “Mouse genome database ID” field under “Curated by the HGNC”
https://www.genenames.org/download/custom/
2.)	RGD (human-rat): use custom download page at HGNC and choose “Rat genome database ID (supplied by RGD)’ field under “Downloaded from external sources”
https://www.genenames.org/download/custom/
3.)	ZFIN: https://zfin.org/downloads
zebrafish-human: https://zfin.org/downloads/human_orthos.txt
zebrafish-mouse: https://zfin.org/downloads/mouse_orthos.txt
zebrafish-fruit fly: https://zfin.org/downloads/fly_orthos.txt

3.	UniProt GeneID mapping files https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/
Next are the example for the species requested for DIOPT assembly by AGR
Human: https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000005640_9606.idmapping
Mouse:
https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000000589_10090.idmapping
Rat:
https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000002494_10116.idmapping
Xenupus:
https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000008143_8364.idmapping
Zebrafish: https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000000437_7955.idmapping
Fruit fly: https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000000803_7227.idmapping
C elegant:
https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000001940_6239.idmapping
Budding yeast: https://ftp.ebi.ac.uk/pub/databases/reference_proteomes/QfO/Eukaryota/UP000002311_559292.idmapping

4.	Gene information file from NCBI
https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
