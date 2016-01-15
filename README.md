### Table of Contents
- **[RapSi](#rapsi)**  
  - **[Prerequisites](#prerequisites)**  
  - **[Running the program](#running-the-program)**  
    - **[Parameters](#parameters)**  
    - **[Test run](#test-run)**  
  - **[FAQ](#faq)**  
  - **[Citation](#citation)**  

# RapSi
**RapSi** (**Rap**id **Si**milarity search) finds similar proteins 

The search consists of two steps:
- **identification of likely targets** by k-mer search against pre-hashed targets (stored in MySQL or SQLite)
- **alignement** of query and selected targets (BLAT)

RapSi is:
- **fast** & **lightweight**, with multi-core support and memory-optimised, reusing existing sequence tables
- **flexible** toward many DB configurations
- **easy-to-implement** with your existing DB schema
- **straightforward** to integrate with front-end of your service (HTML-formatted output)

## Prerequisites
- Python 2.7+ & Biopython 1.6+ `sudo easy_install -U biopython`
- [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)


## Running the program

### Parameters

#### fasta2hash.py


- Genral options:

- SQL options:

- MySQL options:


#### fasta2hits.py

- Genral options:

- SQL options:

- MySQL options:


### Test run

```bash
# get 100 random sequences
src/fasta2hits.py --random 100 --seqcmd "select protid, seq from protid2seq where protid in (%s)" > test.fa

# search agains db
src/fasta2hits.py -i test.fa -o test.fa.out --seqcmd "select protid, seq from protid2seq where protid in (%s)" 

```

## Citation
Leszek P. Pryszcz and Toni Gabaldón (In preparation) RapSi: rapid and highly flexible similarity search for proteins
