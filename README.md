### Table of Contents
- **[RapSi](#rapsi)**  
  - **[Prerequisites](#prerequisites)**  
  - **[Running the program](#running-the-program)**  
    - **[Parameters](#parameters)**  
    - **[Test run](#test-run)**  
  - **[FAQ](#faq)**  
  - **[Citation](#citation)**  

# RapSi
**RapSi** (**Rap**id **Si**milarity search) finds similar sequences nearly instantly. Major bottleneck of existing tools comes from parsing (slow) or loading database of targets (high memory usage) prior to alignment step. To overcome this bottleneck, we hash target sequences into SQL db. This guarantees nearly instant selection of promising targets, while securing low memory footprint. Only pre-selected targets are later aligned with a query. 

RapSi is:
- **fast** & **lightweight**, with multi-core support and memory-optimised, reusing existing sequence tables, with BGZIP compression support
- **flexible** toward many configurations, sequences can be stored in FASTA file or SQL database
- **easy-to-implement** with your existing DB schema
- **straightforward** to integrate with front-end of your service (HTML-formatted output)

RapSi consists of two main modules:
* **fasta2hash.py** - hash database of targets and store it in MySQL or SQLite
* **fasta2hits.py** - similarity search tool 

In order to use RapSi, first you need to hash your target sequences (`fasta2hash.py`). This takes some time, but it's done once per each set of targets. 

The search consists of two steps (`fasta2hits.py`):
- **identification of likely targets** by k-mer search against pre-hashed targets (stored in MySQL or SQLite)
- **alignement** of query and selected targets (BLAT)

<img src="/docs/rapsi_flowchart.png" width="600">

## Prerequisites
- Python 2.7+ & Biopython 1.6+ `sudo easy_install -U biopython`
- [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)
- bgzip (for running test set)

## Running the program
RapSi input consists of either FASTA-formatted file(s) or MySQL DB connection. First, the program will hash the set of targets (`fasta2hash.py`) and store it in SQLite or MySQL DB (depending on input used).
Later, the target sequences can be search for similarity with FASTA or FASTQ-formatted input files (`fasta2hits.py`). 

### Parameters

#### fasta2hash.py
`fasta2hash.py` hash the database of targets and stores hash information in MySQL/SQLite. By default, words of five amino acids (k<sub>H</sub>=5) sliding by one (s<sub>H</sub>=1) are used. Such parameters provide good equilibrium between sensitivity and speed, because 3.2x10<sup>6</sup> possible words (20<sup>5</sup>) can be easily stored and quickly accessed. For the sake of performance, we ignore too common words ie. found in more than 0.05% target sequences, as these bring little information and may negatively affect performance. 

- Genral options:
```
  -h, --help            show this help message and exit
  -v                    verbose
  --version             show program's version number and exit
  -i [INPUT [INPUT ...]], --input [INPUT [INPUT ...]]
                        fasta file(s) [get seqs from db]
  -k KMER, --kmer KMER  hash length [5]
  -r, --replace         overwrite table if exists
  -s STEP, --step STEP  hash steps [1]
  --kmerfrac KMERFRAC   ignore kmers present in more than [0.05] percent
                        targets
  --dna                 DNA alphabet [amino acids]
  --nprocs NPROCS       no. of threads [4]; NOTE: so far only tempfiles
                        parsing is threaded!
  --tmpdir TMPDIR       temp directory [/tmp]
  --tempfiles TEMPFILES
                        temp files no. [1000]
```

- MySQL/SQLite options:
```
  -d DB, --db DB        database [metaphors_201310]
  -t TABLE, --table TABLE
                        hashtable name [hash2protids]
```

- MySQL options:
```
  -u USER, --user USER  database user [lpryszcz]
  -p PSWD, --pswd PSWD  user password [will prompt if not specified]
  --host HOST           database host [localhost]
  -P PORT, --port PORT  server port [3306]
  -c CMD, --cmd CMD     cmd to get protid, seq from db [SELECT protid, seq
                        FROM protid2seq]
  --dtype DTYPE         numpy data type for protids [uint32] ie. S8 for
                        VARCHAR(8) or uint16 for SMALLINT UNSIGNED
  --notempfile          direct upload (without temp file); NOTE: this may
                        cause MySQL time-out on some servers
```


#### fasta2hits.py
`fasta2hits.py` splits query into set of k-mers, loop-up the database of hashed targets, selects candidates for alignment, perform alignments of selected targets with query sequence and return results. 

- Genral options:
```
  -h, --help            show this help message and exit
  -v                    verbose
  --version             show program's version number and exit
  -i INPUT, --input INPUT
                        fasta stream    [stdin]
  -b BATCH, --batch BATCH
                        query batch     [1]
  --seqformat {fasta,fastq,gb,genbank,embl}
                        input format    [auto]
  -X, --rapsiX          6-frames translation of nucleotide query
  --dna                 DNA alphabet    [amino acids]
  -o OUTPUT, --output OUTPUT
                        output stream   [stdout]
  --html                return HTML     [txt]
  --no_query            no query in output
  --link LINK           add html link matches []
  -l LIMIT, --limit LIMIT
                        stop after      [all]
  --random RANDOM       return N random sequence(s) and exit [0]
  --seqcmd SEQCMD       SQL to fetch sequences [select protid, seq from protid2seq where p.protid in (%s)]
```

- Similarity search options:
```
  -k KMER, --kmer KMER  hash length     [5]
  --blatpath BLATPATH   BLAT path       [blat]
  --tmpdir TMPDIR       TEMP path       [.]
  -s STEP, --step STEP  hash steps      [5]
  --seqlimit SEQLIMIT   max. seqs to retrieve [100]
  -n [SAMPLING [SAMPLING ...]], --sampling [SAMPLING [SAMPLING ...]]
                        sample n kmers per query [100]
  --ifrac IFRAC         min. identity of best hit [0.3]
  --ofrac OFRAC         min. overlap of best hit [0.3]
```

- MySQL/SQLite options:
```
  -d DB, --db DB        database        [metaphors_201601]
  -t TABLE, --table TABLE
                        hashtable name  [hash2protids]
```

- MySQL options:
```
  -u USER, --user USER  database user   [script]
  -p PSWD, --pswd PSWD  user password   [python]
  -P PORT, --port PORT  server port     [4042]
  --host HOST           database host   [cgenomics.crg.es]
  --dtype DTYPE         numpy data type for protids [uint32] ie. S8 for VARCHAR(8) or uint16 for SMALLINT UNSIGNED
  -z DBLENGTH, --dblength DBLENGTH
                        database length for E-value calculation [0]
```

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
