### Table of Contents
- **[RapSi](#rapsi)**  
  - **[Prerequisites](#prerequisites)**  
  - **[Running the program](#running-the-program)**  
    - **[Parameters](#parameters)**  
    - **[Test run](#test-run)**  
  - **[FAQ](#faq)**  
  - **[Citation](#citation)**  

# RapSi
**RapSi** (**Rap**id **Si**milarity search) finds similar proteins nearly instantly. Major bottleneck of existing tools comes from parsing (slow) or loading database of targets (high memory usage) prior to alignment step. To overcome this bottleneck, we hash target sequences into SQL db. This guarantees nearly instant selection of promising targets, while securing low memory footprint. Only pre-selected targets are later aligned with a query. 

RapSi is:
- **fast** & **lightweight**, with multi-core support and memory-optimised, reusing existing sequence tables, BGZIP compression support
- **flexible** toward many DB configurations
- **easy-to-implement** with your existing DB schema
- **straightforward** to integrate with front-end of your service (HTML-formatted output)

RapSi consists of two main modules:
* fasta2hash.py - hash database of targets and store in MySQL or SQLite
* fasta2hits.py - similarity search tool that 

In order to use RapSi, first you need to hash your target sequences (fasta2hash.py). This takes some time, but it's done once per each set of targets. 

The search consists of two steps (fasta2hits.py):
- **identification of likely targets** by k-mer search against pre-hashed targets (stored in MySQL or SQLite)
- **alignement** of query and selected targets (BLAT)

![Flowchart](/docs/rapsi_flowchart.png)

## Prerequisites
- Python 2.7+ & Biopython 1.6+ `sudo easy_install -U biopython`
- [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)
- bgzip (for running test set)

## Running the program

### Parameters

#### fasta2hash.py


- Genral options:

- SQL options:

- MySQL options:


By default, we use words of five amino acids (k<sub>H</sub>=5) sliding by one (s<sub>H</sub>=1). Such parameters provide good equilibrium between sensitivity and speed, because 3.2x10<sup>6</sup> possible words (20<sup>5</sup>) can be easily stored and quickly accessed. For the sake of performance, we ignore too common words ie. found in more than 0.05% target sequences, as these bring little information and may negatively affect performance. 

#### fasta2hits.py

- Genral options:

- SQL options:

- MySQL options:


### Test run

```bash
# get uniprot database and bzgip compress it
wget -O- ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz | zcat | bgzip > sprot.gz

# hash database [this make take some minutes]
src/fasta2hash.py -v -i test/sprot.gz -d test/sprot.gz.db3

# get 100 random sequences from db
src/fasta2hits.py -v --random 100 -d test/sprot.gz.db3 > test/test.fa

# search agains db
src/fasta2hits.py -v -i test/test.fa -d test/sprot.gz.db3 -o test/test.fa.out

# check if all targets aligned by comparing query and target ID
awk '$1==$2' test/test.fa.out | wc -l
```

## Citation
Leszek P. Pryszcz and Toni Gabaldón (In preparation) RapSi: rapid and highly flexible similarity search for proteins
