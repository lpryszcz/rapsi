#!/usr/bin/env python
desc="""Generate hash table and load it to db.
Note, you will need ~16GB RAM to hash 13M proteins with kmer=5 & step=1.
MySQLdb connect using SSCursor (server side cursor), as default takes a lot
of RAM. Isn't that slower??
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 13/11/2013
"""

import os, sys, time, zlib
import getpass, resource, subprocess
import array, MySQLdb, MySQLdb.cursors
import sqlite3, tempfile
from datetime import datetime
from Bio import SeqIO, bgzf
from Bio.Alphabet import generic_protein

def adapt_tuple(xlist):
    return array.array('I', [int(x) for x in xlist]).tostring()
    
sqlite3.register_adapter(tuple, adapt_tuple)
sqlite3.register_converter("tuple", array.array.tolist)

aminos   = 'ACDEFGHIKLMNPQRSTVWY'
aminoset = set(aminos)

ALPHABET = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'

def b62encode(n, alphabet=ALPHABET):
    """Converts a positive integer to a base62 string."""
    if not n:
        return alphabet[0]
    base62 = ""
    while n != 0:
        n, i = divmod(n, len(alphabet))
        base62 = alphabet[i] + base62
    return base62#.zfill(leng)

def b62decode(base62, alphabet=ALPHABET):
    """Convert base62 str to int."""
    n = 0
    for i, char in enumerate(base62, 1):
        n += alphabet.index(char) * (len(alphabet) ** (len(base62) - i))
    return n  

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

def seq2mers(seq, kmer, step, aminoset=set(aminos)):
    """Kmers generator for seq"""
    mers = set()
    #skip first amino; usually M
    for s in xrange(1, len(seq)-kmer, step):
        mer = seq[s:s+kmer]
        #skip mer containing non-standard aminos
        if not aminoset.issuperset(mer):
            continue
        mers.add(mer)
    return mers
        
def insert_offset_data(cur, i, fi, psbyte, pebyte, headerlen):
    """Insert record offset info into table."""
    #with previous start and previous end - length of current header
    ebyte = pebyte - headerlen
    cur.execute("INSERT INTO offset_data VALUES (?, ?, ?, ?)",(i, fi, psbyte, ebyte-psbyte))
    #store ebyte
    return ebyte
    
def fasta_parser(fastas, cur, verbose):
    """Fasta iterator returning i, seqid as base64 and sequence str.
    Index sequence using proprietrary parser.
    """
    if verbose:
        sys.stderr.write("Indexing and hashing sequences...\n")
    #parse fasta
    i = 0
    for fi, fn in enumerate(fastas):
        #add file to db
        cur.execute("INSERT INTO file_data VALUES (?, ?)",(fi, fn))
        #get handle and start byte
        if fn.endswith('.gz'):
            handle = bgzf.open(fn)
        else:
            handle = open(fn)
        psbyte = pebyte = 0
        for i, r in enumerate(SeqIO.parse(handle, 'fasta'), i+1):
            #if i>10000: break
            sbyte = handle.tell()
            #store previous entry
            if pebyte:
                psbyte = insert_offset_data(cur, i-1, fi, psbyte, pebyte, len(r.description)+2)
            pebyte = handle.tell()
            yield i, i, str(r.seq) #b62encode(i)
        #store last entry
        if pebyte:
            insert_offset_data(cur, i, fi, psbyte, pebyte, 0)
    #fill metadata
    cur.executemany("INSERT INTO meta_data VALUES (?, ?)",(('count', i),('format','fasta')))
    cur.connection.commit()

def db_seq_parser(cur, cmd, verbose):
    """Database iterator returning i, seqid and seq"""
    l = cur.execute(cmd) 
    if verbose:
        sys.stderr.write("Hashing sequences...\n")
    #parse seqs
    for i, (seqid, seq) in enumerate(cur, 1):
        yield i, seqid, seq
        
def hash_sequences(parser, kmer, step, seqlimit, verbose):
    """Parse input fasta and generate hash table."""
    #open tempfile for each two first aminos
    keylen = 2
    files = {int2mer(i, aminos, keylen): tempfile.TemporaryFile(dir=os.path.curdir) for i in xrange(len(aminos)**keylen)}
    #hash sequences
    i = 0
    for i, seqid, seq in parser:
        if verbose and not i%10e2:
            sys.stderr.write(" %s\r" % i)
        for mer in seq2mers(seq, kmer, step, aminoset):
            files[mer[:keylen]].write("%s\t%s\n"%(mer, seqid))
    sys.stderr.write(" %s\n" % i)
    return files

def parse_tempfiles(files, seqlimit):
    """Generator of mer, protids from each tempfile"""
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
            else:
                mer2protids[mer] = [protid]
        for mer, protids in mer2protids.iteritems():
            yield mer, protids
    
def upload(files, db, host, user, pswd, table, kmer, seqlimit, verbose):        
    """Load to database"""
    #store to database
    if verbose:
        sys.stderr.write("Uploading to database...\n")
    rows = 0
    #fn = "%s.load.%s.txt" % (db, time.time())
    #out = open(fn, "w")
    discarded = 0
    cnx = MySQLdb.connect(db=db, host=host, user=user, passwd=pswd, local_infile = 1)
    cur = cnx.cursor()
    for mer, protids in parse_tempfiles(files, seqlimit):
        if not protids:
            discarded += 1
            continue
        #out.write("%s\t%s\n"%(mer, " ".join(protids)))
        rows += 1
        cur.execute("INSERT INTO `"+table+"` VALUES (%s, %s)",\
                    (mer, buffer(zlib.compress(" ".join(protids)))))        
    '''#close output before uploading!
    out.close()
    #connect & load data
    cnx = MySQLdb.connect(db=db, host=host, user=user, passwd=pswd, local_infile = 1)
    cur = cnx.cursor()
    cmd = "LOAD DATA LOCAL INFILE '%s' INTO TABLE %s"%(fn, table)
    if verbose:
        sys.stderr.write(" %s\n"%cmd)
    cur.execute(cmd)
    #compress table and add index
    #cmd = "myisampack -s %s"'''
    if verbose:
        sys.stderr.write(" %s rows uploaded!\n"%rows)
    cmd = "ALTER TABLE `%s` ADD KEY `idx_hash` (`hash`)" % table
    if verbose:
        sys.stderr.write("Creating index: %s\n"%cmd)
    cur.execute(cmd)
    return discarded
    
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
    #create table #index added after compression 
    #cmd = 'CREATE TABLE `%s` (`hash` CHAR(%s), `protids` MEDIUMTEXT, KEY `idx_hash` (`hash`)) ENGINE=MyISAM' % (table, kmer)
    cmd = 'CREATE TABLE `%s` (`hash` CHAR(%s), `protids` BLOB) ENGINE=MyISAM' % (table, kmer)
    if verbose:
        sys.stderr.write(" %s\n"%cmd)
    cur.execute(cmd)
    return cur
    
def upload_sqlite(files, cur, table, seqlimit, verbose):
    """Load data from tmpfiles into sqlite3."""
    if verbose:
        sys.stderr.write("Loading into database...\n")
    '''#subprocess for inserting
    args = ["sqlite3", db, ".import /dev/stdin %s"%table]
    if verbose:
        sys.stderr.write(" %s\n" % " ".join(args))
    proc = subprocess.Popen(args, stdin=subprocess.PIPE)
    out  = proc.stdin #'''
    #upload
    discarded = 0
    for mer, protids in parse_tempfiles(files, seqlimit):
        if not protids:
            discarded += 1
            continue 
        #out.write('%s|%s\n' % (mer, " ".join(protids))) '''
        #a = array.array('B', [int(x) for x in protids])
        #cur.execute("INSERT INTO %s VALUES (?,?)"%table,(mer, tuple(protids)))
        cur.execute("INSERT INTO %s VALUES (?,?)"%table,\
                    (mer, buffer(zlib.compress(" ".join(protids)))))
    cur.connection.commit()#'''
    return discarded
    
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
    cnx.text_factory = str #sqlite3.OptimizedUnicode
    cur = cnx.cursor()
    #asyn execute >50x faster ##http://www.sqlite.org/pragma.html#pragma_synchronous 
    cur.execute("PRAGMA synchronous=OFF")
    #prepare tables and indices
    cur.execute("CREATE TABLE meta_data (key TEXT, value TEXT)")
    cur.execute("CREATE TABLE file_data (file_number INTEGER, name TEXT)")
    cur.execute("CREATE TABLE offset_data (key INTEGER PRIMARY KEY, file_number INTEGER, offset INTEGER, length INTEGER)")
    #cur.execute("CREATE UNIQUE INDEX key_index ON offset_data(key)")
    cur.execute("CREATE TABLE %s (hash CHAR(%s) PRIMARY KEY, protids BLOB)" % (table, kmer))
    #cur.execute("CREATE INDEX idx_hash ON %s (hash)"%table)
    return cur  
            
def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", "--input",      nargs="*",
                        help="fasta file(s)   [get seqs from db]")
    parser.add_argument("-k", "--kmer",       default=5, type=int, 
                        help="hash length     [%(default)s]")
    parser.add_argument("-r", "--replace",    default=False, action="store_true",
                        help="overwrite table if exists")
    parser.add_argument("-s", "--step",       default=1, type=int, 
                        help="hash steps      [%(default)s]")
    parser.add_argument("--seqlimit",         default=2000, type=int, 
                        help="ignore too common kmers [%(default)s]")
    sqlopt = parser.add_argument_group('MySQL/SQLite options')
    sqlopt.add_argument("-d", "--db",         default="metaphors_201310", 
                        help="database        [%(default)s]")
    sqlopt.add_argument("-t", "--table",      default='hash2protids4',
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
                        
    o = parser.parse_args()

    #get sequence parser
    if o.input:
        #prepare database
        cur = dbConnect_sqlite(o.db, o.table, o.kmer, o.verbose, o.replace)
        #get parser
        parser    = fasta_parser(o.input, cur, o.verbose)
        #hash seqs
        files     = hash_sequences(parser, o.kmer, o.step, o.seqlimit, o.verbose)
        #upload
        discarded = upload_sqlite(files, cur, o.table, o.seqlimit, o.verbose)
    else:
        #prompt for mysql passwd
        pswd = o.pswd
        if pswd==None:
            pswd = getpass.getpass("Enter MySQL password: ")
        #connect
        cur       = dbConnect(o.db, o.host, o.user, pswd, o.table, o.kmer, o.verbose, o.replace)
        #get parser
        parser    = db_seq_parser(cur, o.cmd, o.verbose)
        #hash seqs
        files     = hash_sequences(parser, o.kmer, o.step, o.seqlimit, o.verbose)
        #upload
        discarded = upload(files, o.db, o.host, o.user, pswd, o.table, o.kmer, o.seqlimit, o.verbose)
    
    sys.stderr.write("hash_discarded: %s memory: %s MB\n" % (discarded, memory_usage())) 
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
