#!/usr/bin/env python
desc="""Generate hash table and load it to db.

CHANGELOG:
v1.3:
- np.array compression (30% reduced hash size)
- batch query
- P & E-value statistics
v1.2:
- DNA support
v1.1:
- implemented nucleotide query (rapsiX)
- FastQ, genbank, embl support
- auto query format recognition
- BGZIP support
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona/Mizerow, 13/11/2013
"""

import os, sys, time
import getpass, resource, subprocess
import MySQLdb, MySQLdb.cursors, sqlite3, tempfile
import numpy as np
from datetime import datetime
from Bio import SeqIO, bgzf

aminos = 'ACDEFGHIKLMNPQRSTVWY'
aminoset = set(aminos)
nucleotides = 'ACGT'
DNAcomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def memory_usage():
    """Return memory usage in MB"""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

def int2mer(meri, aminos, kmer):
    """Return char representation of mer
    0      --> AAAAA
    1      --> AAAAC
    20^5-1 --> YYYYY"""
    mer = []
    while meri / len(aminos):
        mer.append(aminos[meri % len(aminos)])
        meri = meri / len(aminos)
    mer.append(aminos[meri % len(aminos)])
    mer.extend(['A']*(kmer-len(mer)))
    mer.reverse()
    return "".join(mer)

def aaseq2mers(seq, kmer, step, aminoset=set(aminos)):
    """Kmers generator for amino seq"""
    """Kmers generator for seq"""
    return set(seq[s:s+kmer] for s in xrange(1, len(seq)-kmer, step) \
               if aminoset.issuperset(seq[s:s+kmer]))

def reverse_complement(mer):
    """Return DNA reverse complement"""
    return "".join(DNAcomplement[b] for b in reversed(mer))
    
def dnaseq2mers(seq, kmer, step, nucleotideset=set(nucleotides)):
    """Kmers generator for DNA seq"""
    mers = set()
    for s in xrange(0, len(seq)-kmer, step):
        mer = seq[s:s+kmer]
        if not nucleotideset.issuperset(mer):
            continue
        #store reverse complement
        if mer > reverse_complement(mer):
            mer = reverse_complement(mer)
        mers.add(mer)
    return mers

def get_seq_offset_length(handle):
    """Return entry start offset and sequence. BGZIP compatible.
    #http://biopython.org/DIST/docs/api/Bio.File-pysrc.html#_SQLiteManySeqFilesDict.__init__
    """
    soffset = eoffset = handle.tell()
    seq = []
    for line in handle:
        if line.startswith(">"):
            if seq:
                yield "".join(seq), soffset, length
            seq = []
            length = len(line)
            soffset = eoffset
        else:
            seq.append(line.strip())
            length += len(line)
            eoffset = handle.tell()
    if seq:
        yield "".join(seq), soffset, length
    
def fasta_parser(fastas, cur, verbose):
    """Fasta iterator returning i, seqid as base64 and sequence str.
    Index sequence using proprietrary parser.
    Handle bgzip compressed files.
    """
    if verbose:
        sys.stderr.write("[%s] Indexing and hashing sequences...\n"%datetime.ctime(datetime.now()))
    #parse fasta
    i = 0
    seqlen = 0
    for fi, fn in enumerate(fastas):
        #add file to db
        cur.execute("INSERT INTO file_data VALUES (?, ?)",(fi, fn))
        #get handle and start byte
        if fn.endswith('.gz'):
            handle = bgzf.open(fn)
        else:
            handle = open(fn)
        #parse entries
        for i, (seq, offset, elen) in enumerate(get_seq_offset_length(handle), i+1):
            if i>100000: break
            seqlen += len(seq)
            cur.execute("INSERT INTO offset_data VALUES (?, ?, ?, ?)",\
                        (i, fi, offset, elen))
            yield i, i, seq
    if verbose:
        sys.stderr.write(" %s letters in %s sequences\n"%(seqlen, i))
    #fill metadata
    cur.executemany("INSERT INTO meta_data VALUES (?, ?)",\
                    (('count', i),('format','fasta'),('dblength',seqlen)))
    #and commit changes
    cur.connection.commit()

def db_seq_parser(cur, cmd, verbose):
    """Database iterator returning i, seqid and seq"""
    l = cur.execute(cmd) 
    if verbose:
        sys.stderr.write("[%s] Hashing sequences...\n"%datetime.ctime(datetime.now()))
    seqlen = 0
    #parse seqs
    for i, (seqid, seq) in enumerate(cur, 1):
        seqlen += len(seq)
        yield i, seqid, seq
    if verbose:
        sys.stderr.write(" %s letters in %s sequences\n"%(seqlen, i))
        
def hash_sequences(parser, kmer, step, dna, verbose):
    """Parse input fasta and generate hash table."""
    #get alphabet
    if dna:
        alphabet = nucleotides
        seq2mers = dnaseq2mers
        keylen   = 4
    else:
        alphabet = aminos
        seq2mers = aaseq2mers
        keylen   = 2
    alphabetset = set(alphabet)
    #open tempfile for each two first aminos
    files = {} 
    for i in xrange(len(alphabet)**keylen):
        files[int2mer(i, alphabet, keylen)] = tempfile.TemporaryFile(dir=os.path.curdir)  
    #hash sequences
    i = 0
    for i, seqid, seq in parser:
        if verbose and not i%10e2:
            sys.stderr.write(" %s\r" % i)
        for mer in seq2mers(seq, kmer, step, alphabetset):
            files[mer[:keylen]].write("%s\t%s\n"%(mer, seqid))
    return files, i

def parse_tempfiles(files, seqlimit, verbose=1):
    """Generator of mer, protids from each tempfile"""
    i = discarded = 0
    for f in files.itervalues():
        f.flush()
        f.seek(0)
        mer2protids = {}
        for l in f:
            mer, protid = l[:-1].split('\t')
            if mer in mer2protids:
                if mer2protids[mer]:
                    mer2protids[mer].append(protid)
                    if len(mer2protids[mer]) > seqlimit:
                        mer2protids[mer] = None
                        discarded += 1
            else:
                mer2protids[mer] = [protid]
        for i, (mer, protids) in enumerate(mer2protids.iteritems(), i+1):
            if not i%10000:
                sys.stderr.write("  %s \r"%i)
            yield mer, protids
    info = " %s hash uploaded; discarded: %s; memory: %s MB\n"
    sys.stderr.write(info%(i-discarded, discarded, memory_usage())) 
            
def upload(files, db, host, user, pswd, table, seqlimit, dtype, verbose, \
           sep = "..|..", end = "..|.\n"):
    """Load to database."""
    args = ["mysql", "-vvv", "-h", host, "-u", user, db, "-e", \
            "LOAD DATA LOCAL INFILE '/dev/stdin' INTO TABLE `%s` FIELDS TERMINATED BY '%s' LINES TERMINATED BY '%s'"%(table, sep, end)]
    if verbose:
        info = "[%s] Uploading to database...\n %s\n" 
        sys.stderr.write(info%(datetime.ctime(datetime.now()), " ".join(args)))
    if pswd:
        args.append("-p%s"%pswd)
    proc = subprocess.Popen(args, stdin=subprocess.PIPE, stderr=sys.stderr)
    out  = proc.stdin
    i = discarded = 0
    for i, (mer, protids) in enumerate(parse_tempfiles(files, seqlimit), 1):
        if not protids:
            discarded += 1
            continue
        out.write("%s%s%s%s"%(mer, sep, np.array(protids, dtype=dtype).tostring(), end))
        if verbose and not i%10000:
            sys.stderr.write("  %s processed; %s discarded\r" % (i, discarded))
    return i, discarded
        
def batch_insert(files, cur, table, seqlimit, verbose, dtype="uint32", rep="?"):
    """Insert to database."""
    if verbose:
        sys.stderr.write("[%s] Uploading...\n"%datetime.ctime(datetime.now()))
    cmd = 'INSERT INTO %s VALUES (%s, %s)'%(table, rep, rep)
    cur.executemany(cmd, ((mer, buffer(np.array(protids, dtype=dtype).tostring())) \
                    for mer, protids in parse_tempfiles(files, seqlimit) if protids))
    #commit
    cur.execute("CREATE INDEX `idx_hash` ON `%s` (hash)"%table)
    cur.connection.commit()
    
def dbConnect(db, host, user, pswd, table, kmer, verbose, replace):
    """Get connection and create empty table.
    Exit if table exists and no overwritting."""
    if verbose:
        sys.stderr.write("Connecting to %s as %s...\n" %(db, user))
    cnx = MySQLdb.connect(db=db, host=host, user=user, passwd=pswd, \
                          cursorclass=MySQLdb.cursors.SSCursor)
    cur = cnx.cursor()
    #check if table exists
    cur.execute('SHOW TABLES')
    tables = set(t for t, in cur.fetchall())
    if table in tables:
        if replace:
            cmd = 'DROP TABLE `%s`'%table
            if verbose:
                sys.stderr.write(" %s\n"%cmd)
            cur.execute(cmd)
        else:
            sys.exit('Table %s already exists!'%table)
    #create table #index added after compression PRIMARY KEY
    cmd = 'CREATE TABLE `%s` (`hash` CHAR(%s), `protids` BLOB) ENGINE=MyISAM'%(table, kmer)
    if verbose:
        sys.stderr.write(" %s\n"%cmd)
    cur.execute(cmd)
    #cur.execute("CREATE INDEX `idx_hash` ON `%s` (hash)"%table)
    return cur

def dbConnect_sqlite(db, table, kmer, verbose, replace):
    """Get connection and create empty table.
    Exit if table exists and no overwritting."""
    if verbose:
        sys.stderr.write("Preparing database: %s...\n" %db)
    #check if db exists
    if os.path.isfile(db):
        #replace or exit
        if replace:
            os.unlink(db)            
        else:
            sys.exit(" Database %s already exists!"%db)
    #connect
    cnx = sqlite3.connect(db)
    #enable utf8 handling
    cnx.text_factory = str 
    cur = cnx.cursor()
    ##bullshit? no time diff
    #asyn execute >50x faster ##http://www.sqlite.org/pragma.html#pragma_synchronous 
    #http://stackoverflow.com/questions/1711631/how-do-i-improve-the-performance-of-sqlite
    cur.execute("PRAGMA synchronous=OFF")
    cur.execute("PRAGMA journal_mode = MEMORY")
    #prepare tables and indices PRIMARY KEY
    cur.execute("CREATE TABLE meta_data (key TEXT, value TEXT)")
    cur.execute("CREATE TABLE file_data (file_number INTEGER, name TEXT)")
    cur.execute("CREATE TABLE offset_data (key INTEGER PRIMARY KEY, file_number INTEGER, offset INTEGER, length INTEGER)")
    cur.execute("CREATE TABLE %s (hash CHAR(%s), protids array)"%(table, kmer))
    return cur
            
def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.3b')   
    parser.add_argument("-i", "--input",      nargs="*",
                        help="fasta file(s)   [get seqs from db]")
    parser.add_argument("-k", "--kmer",       default=5, type=int, 
                        help="hash length     [%(default)s]")
    parser.add_argument("-r", "--replace",    default=False, action="store_true",
                        help="overwrite table if exists")
    parser.add_argument("-s", "--step",       default=1, type=int, 
                        help="hash steps      [%(default)s]")
    parser.add_argument("--kmerfrac",         default=0.01, type=float, 
                        help="ignore kmers present in more than [%(default)s] percent targets")
    parser.add_argument("--dna",              default=False, action='store_true',
                        help="DNA alphabet    [amino acids]")
    sqlopt = parser.add_argument_group('MySQL/SQLite options')
    sqlopt.add_argument("-d", "--db",         default="metaphors_201310", 
                        help="database        [%(default)s]")
    sqlopt.add_argument("-t", "--table",      default='hash2protids',
                        help="hashtable name  [%(default)s]")
    mysqlo = parser.add_argument_group('MySQL options')
    mysqlo.add_argument("-u", "--user",       default="lpryszcz", 
                        help="database user   [%(default)s]")
    mysqlo.add_argument("-p", "--pswd",       default=None, 
                        help="user password   [will prompt if not specified]")
    mysqlo.add_argument("--host",             default='localhost', 
                        help="database host   [%(default)s]")
    mysqlo.add_argument("-c", "--cmd",        default="SELECT protid, seq FROM protid2seq", 
                        help="cmd to get protid, seq from db [%(default)s]")
    mysqlo.add_argument("--dtype",            default="uint32", 
                        help="numpy data type for protids [%(default)s] ie. S8 for VARCHAR(8) or uint16 for SMALLINT UNSIGNED")
                        
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    #get sequence parser
    if o.input:
        #prepare database
        cur = dbConnect_sqlite(o.db, o.table, o.kmer, o.verbose, o.replace)
        #get parser
        parser = fasta_parser(o.input, cur, o.verbose)
        #hash seqs
        files, targets = hash_sequences(parser, o.kmer, o.step, o.dna, o.verbose)
        seqlimit = o.kmerfrac * targets / 100
        #upload
        batch_insert(files, cur, o.table, seqlimit, o.verbose)
    else:
        #prompt for mysql passwd
        pswd = o.pswd
        if pswd==None:
            pswd = getpass.getpass("Enter MySQL password: ")
        #connect
        cur = dbConnect(o.db, o.host, o.user, pswd, o.table, o.kmer, o.verbose, \
                        o.replace)
        #get parser
        parser = db_seq_parser(cur, o.cmd, o.verbose)
        #hash seqs
        files, targets = hash_sequences(parser, o.kmer, o.step, o.dna, o.verbose)
        seqlimit = o.kmerfrac * targets / 100
        #upload
        upload(files, o.db, o.host, o.user, pswd, o.table, seqlimit, o.dtype, o.verbose)
        #batch_insert(files, cur, o.table, seqlimit, o.verbose, o.dtype, "%s")
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
