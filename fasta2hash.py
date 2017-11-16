#!/usr/bin/env python
desc="""Generate hash table and load it to db.

Using reduced amino alphabet and recoding DNA bases as dinucleotides.
In addition, sequences longer than 40K are divided into chunks.

Limits:
- uint32 seq chunks
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona/Mizerow/Warsaw, 13/11/2013
"""

import os, sys, time
import commands, gc, getpass, resource, subprocess
import MySQLdb, MySQLdb.cursors, sqlite3, tempfile
import numpy as np
from datetime import datetime
from Bio import SeqIO, bgzf
from multiprocessing import Pool, Process, Queue

aminos = 'ACDEFGHIKLMNPQRSTVWY'
aminoset = set(aminos)
nucleotides = 'ACGT'
DNAcomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

# http://www.rpgroup.caltech.edu/publications/supplements/alphabets/HP/Welcome.html
## MFILV, ACW, YQHPGTSN, RK, DE
mj5 = {'A': '1', 'C': '1', 'E': '4', 'D': '4', 'G': '2', 'F': '0', 'I': '0', 'H': '2', 'K': '3', 'M': '0', 'L': '0', 'N': '2', 'Q': '2', 'P': '2', 'S': '2', 'R': '3', 'T': '2', 'W': '1', 'V': '0', 'Y': '2'}
## MFILV, A, C, WYQHPGTSN, RK, DE
mj6 = {'A': '1', 'C': '2', 'E': '5', 'D': '5', 'G': '3', 'F': '0', 'I': '0', 'H': '3', 'K': '4', 'M': '0', 'L': '0', 'N': '3', 'Q': '3', 'P': '3', 'S': '3', 'R': '4', 'T': '3', 'W': '3', 'V': '0', 'Y': '3'}
reduced_alphabet = mj6

# recode di-nucleotides as base16 (0-F) characters
di4 = {'AA': '0', 'AC': '1', 'GT': 'B', 'AG': '2', 'CC': '5', 'CA': '4', 'CG': '6', 'TT': 'F', 'GG': 'A', 'GC': '9', 'AT': '3', 'GA': '8', 'TG': 'E', 'TA': 'C', 'TC': 'D', 'CT': '7'}
reduced_dna_alphabet = di4

def memory_usage():
    """Return memory usage in MB"""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024

def encode(n, alphabet, length=5, base=""):
    """Converts a positive integer to a string based in given alphabet."""
    while n != 0:
        n, i = divmod(n, len(alphabet))
        base = alphabet[i] + base
    return base.rjust(length).replace(' ',alphabet[0])
        
def decode(base, alphabet, n=0):
    """Convert string based in given alphabet to int."""
    for i, char in enumerate(base, 1):
        n += alphabet.index(char) * (len(alphabet) ** (len(base) - i))
    return n  
    
def aaseq2mers(seq, kmer, step, aminoset=set(reduced_alphabet.values())):
    """Kmers generator for amino seq"""
    # reduce sequence
    seq = "".join(reduced_alphabet[aa] if aa in reduced_alphabet else aa for aa in seq)
    mers = set(seq[s:s+kmer] for s in xrange(0, len(seq)-kmer, step))
    # return int representation of mer starting from 1
    return [int(mer, base=len(aminoset))+1 for mer in filter(lambda x: x.isdigit(), mers)]

def aaseq2mers0(seq, kmer, step, aminoset=set(reduced_alphabet.values())):
    """Kmers generator for amino seq"""
    mers = set(seq[s:s+kmer] for s in xrange(0, len(seq)-kmer, step))
    for mer in mers:
        # get mer as int with given base 
        imer = "".join(map(str, (reduced_alphabet[aa] for aa in mer if aa in reduced_alphabet)))
        # catch wrong kmers
        if len(imer)!=len(mer):
            continue
        # avoid 0
        yield int(imer, base=len(aminoset))+1

def reverse_complement(mer):
    """Return DNA reverse complement"""
    return "".join(DNAcomplement[b] for b in reversed(mer) if b in DNAcomplement)
    
def dnaseq2mers(seq, kmer, step, baseset=set(reduced_dna_alphabet.values())):
    """Kmers generator for DNA seq"""
    # di-nucleotides - shift by one to make sure both variants will be present
    mers = set(seq[s:s+kmer] for s in xrange(0, len(seq)-kmer, step))#.union(seq[s:s+kmer] for s in xrange(1, len(seq)-kmer, step))
    for mer in mers:
        #store reverse complement
        if mer > reverse_complement(mer):
            mer = reverse_complement(mer)
        # get mer as int with given base (dinucleotide)
        imer = "".join(map(str, (reduced_dna_alphabet[mer[i:i+2]] for i in range(0, kmer, 2) if mer[i:i+2] in reduced_dna_alphabet)))
        # catch wrong kmers
        if len(imer)*2!=kmer:
            continue
        # avoid 0
        yield int(imer, base=len(baseset))+1

def get_seq_offset_length(handle, chunksize=50000):
    """Return name, sequence, entry start offset and length. BGZIP compatible. 
    #http://biopython.org/DIST/docs/api/Bio.File-pysrc.html#_SQLiteManySeqFilesDict.__init__
    """
    name, pstart, length, seq = '', 0, 0, []
    #for line in handle:
    while True: # buffering = 1 doesn't work, so need while
        line = handle.readline()
        if not line:
            break
        if line.startswith(">"):
            if seq:
                _seq = "".join(seq)
                yield "%s:%s-%s"%(name, pstart, pstart+len(_seq)), _seq, soffset, length
            # new name and reset sequence
            name, pstart, length, seq = line[1:].split()[0], 0, 0, []
            soffset = handle.tell()
        else:  
            seq.append(line.strip())
            length += len(line)
            # report seq chunk
            if sum(map(len, seq))>chunksize:
                _seq = "".join(seq)
                yield "%s:%s-%s"%(name, pstart, pstart+len(_seq)), _seq, soffset, length
                pstart += len(_seq)
                seq, length = [], 0
                soffset = handle.tell()
                
    if seq:
        _seq = "".join(seq)
        yield "%s:%s-%s"%(name, pstart, pstart+len(_seq)), _seq, soffset, length
    
def fasta_parser(fastas, cur, verbose):
    """Fasta iterator returning seqid as base64 and sequence str.
    Index sequence using proprietrary parser.
    Handle bgzip compressed files.
    """
    if verbose:
        sys.stderr.write("[%s] Hashing and indexing sequences...\n"%datetime.ctime(datetime.now()))
    #parse fasta
    i = 0
    seqlen = 0
    cmd = "INSERT INTO offset_data VALUES (?, ?, ?, ?, ?)"
    for fi, fn in enumerate(fastas):
        #add file to db
        cur.execute("INSERT INTO file_data VALUES (?, ?)", (fi, fn))
        #get handle and start byte
        if fn.endswith('.gz'):
            handle = bgzf.open(fn) 
        else:
            handle = open(fn)
        #parse entries
        for i, (name, seq, offset, elen) in enumerate(get_seq_offset_length(handle), i+1):
            #if i>10**6: break
            seqlen += len(seq)
            cur.execute(cmd, (i, fi, name, offset, elen))
            yield i, seq
    if verbose:
        sys.stderr.write(" %s letters in"%seqlen)
    #fill metadata
    cur.executemany("INSERT INTO meta_data VALUES (?, ?)", \
                    (('count', i), ('format', 'fasta'), ('dblength', seqlen)))
    #and commit changes
    cur.connection.commit()

def db_seq_parser(cur, cmd, verbose):
    """Database iterator returning seqid and seq"""
    l = cur.execute(cmd) 
    if verbose:
        sys.stderr.write("[%s] Hashing sequences...\n"%datetime.ctime(datetime.now()))
    seqlen = 0
    #parse seqs
    for seqid, seq in cur:
        seqlen += len(seq)
        yield seqid, seq
    if verbose:
        sys.stderr.write(" %s letters in"%seqlen)

def worker1(inQ, parser, nproc):
    """Reading sequences from parser and feeding them to queue"""
    for seqid, seq in parser:
        inQ.put((seqid, seq))
    for i in range(nproc):
        inQ.put(None)
        
def worker2(inQ, outQ, kmer, step, seq2mers, alphabetset):
    """Hashing sequences from the queue"""
    for data in iter(inQ.get, 1):
        if not data:
            break
        seqid, seq = data
        # get mers from upper-case seq
        outQ.put((seqid, seq2mers(seq.upper(), kmer, step, alphabetset)))
    outQ.put(None)
    
def hash_sequences_multi(parser, kmer, step, dna, kmerfrac, tmpdir='/tmp', \
                         tmpfiles=1000, nproc=4, verbose=1): 
    """Parse input fasta and generate hash table."""
    #get alphabet
    if dna:
        alphabet = reduced_dna_alphabet.values()
        seq2mers = dnaseq2mers
        merspace = len(alphabet)**kmer/2 # reverse complement
    else:
        alphabet = reduced_alphabet.values() # aminos
        seq2mers = aaseq2mers
        merspace = len(alphabet)**kmer
    alphabetset = set(alphabet)
    if verbose:
        info = "[%s] Preparing %s temporary files...\n"
        sys.stderr.write(info%(datetime.ctime(datetime.now()), tmpfiles))
    # open tempfiles
    files = [tempfile.NamedTemporaryFile(dir=tmpdir, delete=0) for i in xrange(tmpfiles)]
    # start workers 
    inQ, outQ = Queue(1000), Queue(1000)
    # 1 thread for seq reading
    Process(target=worker1, args=(inQ, parser, nproc)).start()
    # multiple threads for hashing
    for i in range(nproc):
        Process(target=worker2, args=(inQ, outQ, kmer, step, seq2mers, alphabetset)).start()        
    # hash sequences
    i = 0
    stops = 0
    for data in iter(outQ.get, 1):
        if not data:
            stops += 1
            if stops==nproc:
                break
            continue
        i += 1
        if verbose and not i%1e3:
            sys.stderr.write(" %s\r"%i)
            #break
        seqid, mers = data
        for mer in mers:
            files[mer%len(files)].write("%s\t%s\n"%(mer, seqid))
    # set seqlimit only if >10k sequences
    seqlimit = int(round(kmerfrac * i / 100)) if i > 1e5 else i
    info = " %s sequences / chunks [memory: %s MB]\n Setting seqlimit to: %s\n"
    sys.stderr.write(info%(i, memory_usage(), seqlimit))
    return files, seqlimit

def hash_sequences_single(parser, kmer, step, dna, kmerfrac, tmpdir='/tmp', \
                          tmpfiles=1000, nproc=1, verbose=1):
    """Parse input fasta and generate hash table."""
    #get alphabet
    if dna:
        alphabet = reduced_dna_alphabet.values()
        seq2mers = dnaseq2mers
        merspace = len(alphabet)**kmer/2 # reverse complement
    else:
        alphabet = reduced_alphabet.values() # aminos
        seq2mers = aaseq2mers
        merspace = len(alphabet)**kmer
    alphabetset = set(alphabet)
    if verbose:
        info = "[%s] Preparing %s temporary files...\n"
        sys.stderr.write(info%(datetime.ctime(datetime.now()), tmpfiles))
    # open tempfiles
    files = [tempfile.NamedTemporaryFile(dir=tmpdir, delete=0) for i in xrange(tmpfiles)]
    #hash sequences
    i = 0
    for i, (seqid, seq) in enumerate(parser, 1):
        if verbose and not i%1e3:
            sys.stderr.write(" %s\r"%i)
            #break
        # get mers from upper-case seq
        for mer in seq2mers(seq.upper(), kmer, step, alphabetset):
            files[mer%len(files)].write("%s\t%s\n"%(mer, seqid))
    # set seqlimit only if >10k sequences
    seqlimit = int(round(kmerfrac * i / 100)) if i > 1e5 else i
    info = " %s sequences / chunks [memory: %s MB]\n Setting seqlimit to: %s\n"
    sys.stderr.write(info%(i, memory_usage(), seqlimit))
    return files, seqlimit
    
def worker(inQ, outQ, dtype, seqlimit): 
    """Tempfile parsing worker."""
    #py2.6 path
    gc.disable()
    for fn in iter(inQ.get, None):
        mer2protids = {}
        for l in open(fn):
            mer, protid = l[:-1].split('\t')
            if mer in mer2protids:
                if mer2protids[mer]:
                    mer2protids[mer].append(protid)
                    if len(mer2protids[mer])>seqlimit:
                        mer2protids[mer] = []
            else:
                mer2protids[mer] = [protid]
        #convert to np.array compressed string representation
        for mer, protids in mer2protids.iteritems():
            outQ.put((mer, np.array(protids, dtype=dtype).tostring()))
        #remove tmp file
        os.unlink(fn)
        #garbage collection
        del mer2protids
        gc.collect()
    outQ.put(None)

def parse_tempfiles(files, seqlimit, dtype, nprocs=4, verbose=1):
    """Generator of mer, protids from each tempfile"""
    inQ, outQ = Queue(), Queue(100000)
    #populate inQ
    for f in files:
        f.close()
        inQ.put(f.name)
    #start workers
    for i in range(nprocs):
        Process(target=worker, args=(inQ, outQ, dtype, seqlimit)).start()
        inQ.put(None)
    #process pool output
    i = discarded = stop = 0
    for i, data in enumerate(iter(outQ.get, 1), 1-nprocs):
        if not data:
            stop += 1
            if stop==nprocs:
                break
            continue
        (mer, protids) = data
        if not i%10000:
            sys.stderr.write("  %s [memory: %s MB]   \r"%(i, memory_usage()))
        if not protids:
            discarded += 1
            continue
        yield mer, buffer(protids)
    info = " %s hash uploaded; discarded: %s [memory: %s MB]\n"
    sys.stderr.write(info%(i-discarded, discarded, memory_usage())) 
                
def upload(files, db, host, port, user, pswd, table, seqlimit, dtype, nprocs, \
           notempfile=0, tmpdir="./", verbose=1, step=1000):
    """Load to database, optionally through tempfile."""
    cur = _connect(db, host, port, user, pswd)
    cmd = "INSERT INTO "+table+" (hash, protids) VALUES (%s, %s)"
    # execute in batches, to avoid error due to 'max_allowed_packet'
    data = []
    for (mer, protids) in parse_tempfiles(files, seqlimit, dtype, nprocs, verbose):
        data.append((mer, protids))
        if len(data)>step:
            cur.executemany(cmd, data)
            data = []
    if data:
        cur.executemany(cmd, data)
                           
def batch_insert(files, cur, table, seqlimit, nprocs, verbose, \
                 dtype="uint32", rep="?"):
    """Insert to database."""
    if verbose:
        sys.stderr.write("[%s] Uploading...\n"%datetime.ctime(datetime.now()))
    cmd = 'INSERT INTO %s VALUES (%s, %s)'%(table, rep, rep)
    cur.executemany(cmd, parse_tempfiles(files, seqlimit, dtype, nprocs, verbose))
    #commit
    cur.execute("CREATE INDEX `idx_hash` ON `%s` (hash)"%table)
    cur.connection.commit()

def _connect(db, host, port, user, pswd):
    """Return db cursor"""
    cnx = MySQLdb.connect(db=db, host=host, port=port, user=user, passwd=pswd, \
                          cursorclass=MySQLdb.cursors.SSCursor)
    cur = cnx.cursor()
    return cur
    
def dbConnect(db, host, port, user, pswd, table, kmer, verbose, replace):
    """Get connection and create empty table.
    Exit if table exists and no overwritting."""
    if verbose:
        info = "[%s] Connecting to %s @ %s as %s ...\n"
        sys.stderr.write(info%(datetime.ctime(datetime.now()), db, host, user))
    cur = _connect(db, host, port, user, pswd)
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
    # create table #index added after compression PRIMARY KEY
    #cmd = 'CREATE TABLE `%s` (`hash` CHAR(%s), `protids` BLOB, INDEX `idx_hash` (hash)) ENGINE=MyISAM'%(table, kmer)
    ## INT - optimise it
    cmd = 'CREATE TABLE `%s` (`hash` INT(%s), `protids` BLOB, INDEX `idx_hash` (hash)) ENGINE=MyISAM'%(table, kmer)
    if verbose:
        sys.stderr.write(" %s\n"%cmd)
    cur.execute(cmd)
    return cur

def dbConnect_sqlite(db, table, kmer, verbose, replace):
    """Get connection and create empty table.
    Exit if table exists and no overwritting."""
    if verbose:
        info = "[%s] Preparing sqlite3 database: %s ...\n"
        sys.stderr.write(info%(datetime.ctime(datetime.now()), db))
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
    cur.execute("CREATE TABLE offset_data (key INTEGER PRIMARY KEY, file_number INTEGER, seqname TEXT, offset INTEGER, length INTEGER)")
    cur.execute("CREATE TABLE %s (hash INT(%s), protids array)"%(table, kmer))
    return cur
            
def main(): 
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='%(prog)s 1.5a')   
    parser.add_argument("-i", "--input",      nargs="*",
                        help="fasta file(s)   [get seqs from db]")
    parser.add_argument("-k", "--kmer",       default=10, type=int, 
                        help="hash length     [%(default)s]")
    parser.add_argument("-r", "--replace",    default=False, action="store_true",
                        help="overwrite table if exists")
    parser.add_argument("-s", "--step",       default=1, type=int, 
                        help="hash steps      [%(default)s]")
    parser.add_argument("--kmerfrac",         default=0.05, type=float, 
                        help="ignore kmers present in more than [%(default)s] percent targets")
    parser.add_argument("--dna",              default=False, action='store_true',
                        help="DNA alphabet    [amino acids]")
    parser.add_argument("--nprocs",           default=4, type=int, 
                        help="no. of threads  [%(default)s]")
    parser.add_argument("--tmpdir",           default="/tmp",
                        help="temp directory  [%(default)s]")
    parser.add_argument("--tempfiles",        default=1000, type=int, 
                        help="no. of tmp files  [%(default)s]")
    sqlopt = parser.add_argument_group('MySQL/SQLite options')
    sqlopt.add_argument("-d", "--db",         default="", 
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
    mysqlo.add_argument("-P", "--port",       default=4042, type=int, 
                        help="server port     [%(default)s]")
    mysqlo.add_argument("-c", "--cmd",        default="SELECT protid, seq FROM protid2seq", 
                        help="cmd to get protid, seq from db [%(default)s]")
    mysqlo.add_argument("--dtype",            default="uint32", 
                        help="numpy data type for protids [%(default)s] ie. S8 for VARCHAR(8) or uint16 for SMALLINT UNSIGNED")
    mysqlo.add_argument("--notempfile",       default=False, action="store_true", 
                        help="direct upload (without temp file); NOTE: this may cause MySQL time-out on some servers")
                        
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    hash_sequences = hash_sequences_multi
    if o.nprocs<2:
        hash_sequences = hash_sequences_single
        
    # get sequence parser
    if o.input:
        # prepare database
        cur = dbConnect_sqlite(o.db, o.table, o.kmer, o.verbose, o.replace)
        # get parser
        parser = fasta_parser(o.input, cur, o.verbose)
        # hash seqs
        t0 = datetime.now()
        files, seqlimit = hash_sequences(parser, o.kmer, o.step, o.dna, o.kmerfrac, \
                                         o.tmpdir, o.tempfiles, o.nprocs, o.verbose)
        t1 = datetime.now()
        print t1-t0#; return
        # upload
        batch_insert(files, cur, o.table, seqlimit, o.nprocs, o.verbose)
        t2 = datetime.now()
        print t2-t1
    else:
        # prompt for mysql passwd
        pswd = o.pswd
        if pswd==None:
            pswd = getpass.getpass("Enter MySQL password: ")
        # connect
        cur = dbConnect(o.db, o.host, o.port, o.user, pswd, o.table, o.kmer, o.verbose, \
                        o.replace)
        # get parser
        parser = db_seq_parser(cur, o.cmd, o.verbose)
        # hash seqs
        files, seqlimit = hash_sequences(parser, o.kmer, o.step, o.dna, o.kmerfrac, \
                                         o.tmpdir, o.tempfiles, o.nprocs, o.verbose)
        # upload
        upload(files, o.db, o.host, o.port, o.user, pswd, o.table, seqlimit, o.dtype, \
               o.nprocs, o.notempfile, o.tmpdir, o.verbose)
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
