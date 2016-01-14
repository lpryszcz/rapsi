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
                   
def aaseq2mers(seq, kmer, step, aminoset=set(aminos)):
    """Kmers generator for amino seq"""
    return set(seq[s:s+kmer] for s in xrange(0, len(seq)-kmer, step))

def reverse_complement(mer):
    """Return DNA reverse complement"""
    return "".join(DNAcomplement[b] for b in reversed(mer))
    
def dnaseq2mers(seq, kmer, step, nucleotideset=set(nucleotides)):
    """Kmers generator for DNA seq"""
    mers = set()
    for s in xrange(0, len(seq)-kmer, step):
        mer = seq[s:s+kmer]
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
        sys.stderr.write("[%s] Hashing and indexing sequences...\n"%datetime.ctime(datetime.now()))
    #parse fasta
    i = 0
    seqlen = 0
    cmd = "INSERT INTO offset_data VALUES (?, ?, ?, ?)"
    for fi, fn in enumerate(fastas):
        #add file to db
        cur.execute("INSERT INTO file_data VALUES (?, ?)", (fi, fn))
        #get handle and start byte
        if fn.endswith('.gz'):
            handle = bgzf.open(fn)
        else:
            handle = open(fn)
        #parse entries
        for i, (seq, offset, elen) in enumerate(get_seq_offset_length(handle), i+1):
            #if i>10**6: break
            seqlen += len(seq)
            cur.execute(cmd, (i, fi, offset, elen))
            yield i, seq
    if verbose:
        sys.stderr.write(" %s letters in"%seqlen)
    #fill metadata
    cur.executemany("INSERT INTO meta_data VALUES (?, ?)", \
                    (('count', i), ('format', 'fasta'), ('dblength', seqlen)))
    #and commit changes
    cur.connection.commit()

def db_seq_parser(cur, cmd, verbose):
    """Database iterator returning i, seqid and seq"""
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
        
def hash_sequences(parser, kmer, step, dna, kmerfrac, tmpdir='/tmp', \
                   tmpfiles=1000, verbose=1):
    """Parse input fasta and generate hash table."""
    #get alphabet
    if dna:
        alphabet = nucleotides
        seq2mers = dnaseq2mers
        merspace = len(alphabet)**kmer/2 #reverse complement
    else:
        alphabet = aminos
        seq2mers = aaseq2mers
        merspace = len(alphabet)**kmer
    alphabetset = set(alphabet)
    if verbose:
        info = "[%s] Preparing %s temporary files...\n"
        sys.stderr.write(info%(datetime.ctime(datetime.now()), tmpfiles))
    #open tempfiles
    files = [tempfile.NamedTemporaryFile(dir=tmpdir, delete=0) for i in xrange(tmpfiles)]
    #and link every possible kmer to file - should pay off 
    mer2file = {} #py2.6 compatible
    for i in xrange(merspace):
        mer2file[encode(i, alphabet, 5)] = files[i%len(alphabet)]
    #hash sequences
    i = 0
    for i, (seqid, seq) in enumerate(parser, 1):
        if verbose and not i%1e3:
            sys.stderr.write(" %s\r"%i)
        for mer in seq2mers(seq, kmer, step, alphabetset):
            #catch wrong kmers
            if mer in mer2file:
                mer2file[mer].write("%s\t%s\n"%(mer, seqid))
    #get seqlimit
    seqlimit = int(kmerfrac * i / 100)
    info = " %s sequences [memory: %s MB]\n Setting seqlimit to: %s\n"
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
           notempfile=0, tmpdir="./", verbose=1, sep = "....|....", end = "....|....\n"):
    """Load to database, optionally through tempfile."""
    args = ["mysql", "-vvv", "-h", host, "-P", port, "-u", user, db, "-e", \
            "LOAD DATA LOCAL INFILE '/dev/stdin' INTO TABLE `%s` FIELDS TERMINATED BY %s LINES TERMINATED BY %s"%(table, repr(sep), repr(end))]
    args = map(str, args)
    if pswd:
        args.append("-p%s"%pswd)
    #write to mysql directly
    if notempfile:
        p = subprocess.Popen(args, stdin=subprocess.PIPE, stdout=subprocess.PIPE, \
                             stderr=subprocess.PIPE)
        out = p.stdin
        if verbose:
            info = "[%s] Uploading directly to database...\n %s\n"
            sys.stderr.write(info%(datetime.ctime(datetime.now()), \
                                   " ".join(filter(lambda x: not x.startswith("-p"), args))))
    #or through tempfile
    else:
        out = tempfile.NamedTemporaryFile(dir=tmpdir, delete=0)
        if verbose:
            info = "[%s] Storing to tempfile: %s ...\n"
            sys.stderr.write(info%(datetime.ctime(datetime.now()), out.name))
    #write to out
    try:
        out.write("".join("%s%s%s%s"%(mer, sep, protids, end) for mer, protids in \
                          parse_tempfiles(files, seqlimit, dtype, nprocs, verbose)))
        #close out
        out.close()
    #with exception catching
    except Exception, e:
        sys.stderr.write("[ERROR] %s\n"%str(e))
    #start subprocess uploading the data
    if not notempfile:
        #upload from tempfile, instead stdin
        args[10] = args[10].replace('/dev/stdin', out.name)
        if verbose:
            info = "[%s] Uploading to database...\n %s\n"
            sys.stderr.write(info%(datetime.ctime(datetime.now()), \
                                   " ".join(filter(lambda x: not x.startswith("-p"), args))))
        #upload tmpfile
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #wait to finish & check return code
    p.wait()
    if p.returncode:
        info = "[WARNING] Likely error on data upload:\n%s\n"
        sys.stderr.write(info%"\n".join(p.stderr.readlines()))
        if not notempfile:
            sys.stderr.write("NOTE: You can reuse the temp file: %s !\n"%out.name)
    elif not notempfile:
        os.unlink(out.name)
                          
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
    
def dbConnect(db, host, port, user, pswd, table, kmer, verbose, replace):
    """Get connection and create empty table.
    Exit if table exists and no overwritting."""
    if verbose:
        info = "[%s] Connecting to %s @ %s as %s ...\n"
        sys.stderr.write(info%(datetime.ctime(datetime.now()), db, host, user))
    cnx = MySQLdb.connect(db=db, host=host, port=port, user=user, passwd=pswd, \
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
    cmd = 'CREATE TABLE `%s` (`hash` CHAR(%s), `protids` BLOB, INDEX `idx_hash` (hash)) ENGINE=MyISAM'%(table, kmer)
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
    parser.add_argument("--kmerfrac",         default=0.05, type=float, 
                        help="ignore kmers present in more than [%(default)s] percent targets")
    parser.add_argument("--dna",              default=False, action='store_true',
                        help="DNA alphabet    [amino acids]")
    parser.add_argument("--nprocs",           default=4, type=int, 
                        help="no. of threads  [%(default)s]; NOTE: so far only tempfiles parsing is threaded!")
    parser.add_argument("--tmpdir",           default="/tmp",
                        help="temp directory  [%(default)s]")
    parser.add_argument("--tempfiles",        default=1000, type=int, 
                        help="temp files no.  [%(default)s]")
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
    mysqlo.add_argument("-P", "--port",       default=3306, type=int, 
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

    #get sequence parser
    if o.input:
        #prepare database
        cur = dbConnect_sqlite(o.db, o.table, o.kmer, o.verbose, o.replace)
        #get parser
        parser = fasta_parser(o.input, cur, o.verbose)
        #hash seqs
        files, seqlimit = hash_sequences(parser, o.kmer, o.step, o.dna, o.kmerfrac, \
                                         o.tmpdir, o.tempfiles, o.verbose)
        #upload
        batch_insert(files, cur, o.table, seqlimit, o.nprocs, o.verbose)
    else:
        #prompt for mysql passwd
        pswd = o.pswd
        if pswd==None:
            pswd = getpass.getpass("Enter MySQL password: ")
        #connect
        cur = dbConnect(o.db, o.host, o.port, o.user, pswd, o.table, o.kmer, o.verbose, \
                        o.replace)
        #get parser
        parser = db_seq_parser(cur, o.cmd, o.verbose)
        #hash seqs
        files, seqlimit = hash_sequences(parser, o.kmer, o.step, o.dna, o.kmerfrac, \
                                         o.tmpdir, o.tempfiles, o.verbose)
        #upload
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
