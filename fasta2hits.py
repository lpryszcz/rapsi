#!/usr/bin/env python
desc="""Scan hash table and report matches. For searches for relatively similar
sequences, sequence sampling is recommended.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona/Mizerow, 13/11/2013
"""

import commands, gzip, math, os, sys, time
import MySQLdb, random, resource, sqlite3, tempfile
import numpy as np
from datetime import datetime
from Bio import SeqIO, bgzf
from htmlTable import htmlTable
from fasta2hash import dnaseq2mers, aaseq2mers
        
def sqlite2seq(cur, db, protids):
    """Return target fastas for protids from sqlite3."""
    #open target files
    cur.execute("SELECT name FROM file_data")
    files = {} #name: open(os.path.join(os.path.dirname(db), name)) for name, in cur.fetchall()}
    for name, in cur.fetchall():
        # if db is in another directory
        if os.path.isfile(name):
            fpath = name
        else:
            fpath = os.path.join(os.path.dirname(db), name)
        # open fasta file
        if name.endswith('.gz'):
            files[name] = bgzf.open(fpath)
        else:
            files[name] =      open(fpath)

    #get targets
    cmd = """SELECT f.name, offset, length FROM offset_data o JOIN file_data f
    ON o.file_number=f.file_number WHERE key IN (%s)""" % ",".join(str(p) for p in protids)
    cur.execute(cmd)
    targets = []
    for name, offset, length in sorted(cur.fetchall()):
        try:
            files[name].seek(offset)
            targets.append(files[name].read(length))
        except:
            #bgzip sometimes doesn't work at first seek
            sys.stderr.write("[Warning] Cannot fetch sequence for %s at %s + %s bytes\n"%(name, offset, length))
    return "".join(targets)

def mysql2seq(cur, protids, seqcmd):
    """Return target fasta for protids from mysql."""
    protidstext = "'" + "','".join(str(p) for p in protids) + "'"
    cur.execute(seqcmd%protidstext)
    return "".join(">%s\n%s\n"%tup for tup in cur.fetchall())
    
def seq2matches(cur, db, table, seqcmd, qid, qseqs, kmer, step, seqlimit, samplings, \
                seq2mers, dtype, verbose):
    """Return matching protids and sequences"""
    if verbose:
        info = "Parsing %s aminos from %s sequences...\n"
        sys.stderr.write(info % (len("".join(qseqs)), len(qseqs)))
    kmers = [] #set()
    for qseq in qseqs:
        kmers.append(seq2mers(qseq, kmer, step))
    if verbose:
        sys.stderr.write("  %s mers\n"%sum(len(x) for x in kmers))
    fprotids = []
    for sampling in samplings:
        mers = []
        for qmers in kmers:
            mers += list(qmers)[:sampling]
        mers = set(mers)
        if verbose:
            sys.stderr.write("  Sampled %s mers: %s...\n"% (len(mers),", ".join(list(mers)[:6])))
        #get protids for set of mers
        cmd = "select protids from `%s` where hash in ('%s')" % (table, "','".join(mers))    
        cur.execute(cmd)
        protids = {}
        for merprotids, in cur:
            #unpack numpy object
            for protid in np.fromstring(merprotids, dtype=dtype):
                if protid not in protids:
                    protids[protid]  = 1
                else:
                    protids[protid] += 1
        #select upto 100 best seqs with most kmers matching
        fprotids = protids
        if len(fprotids)>seqlimit:
            counts = [0]*(max(protids.itervalues())+1)
            for p, c in protids.iteritems():
                counts[c] += 1
            maxc = 0
            for i, c in enumerate(reversed(counts), 1):
                maxc += c
                if maxc>seqlimit:
                    break
            n = len(counts)-i
            fprotids = filter(lambda x: protids[x]>n, protids)
            if verbose:
                sys.stderr.write("   %8i protids having %s+ kmer matches.\n"% (len(fprotids), n))
        if fprotids:
            break
    #return None if no kmer matches
    if not fprotids:
        return 
    #get sequences
    if verbose:
        sys.stderr.write(" Fetching %s sequences...\n" % len(fprotids))
    #select protids as text or str
    if not os.path.isfile(db): #seqcmd:
        targets = mysql2seq(cur, fprotids, seqcmd)
    else:
        targets = sqlite2seq(cur, db, fprotids)
    return targets
    
def hits2algs(qids, qseqs, matches, blatpath, tmpdir, link, dblength, dna, verbose, \
              Lambda=0.318, K=0.13):
    """Align query with hits and return global algs."""
    if verbose:
        sys.stderr.write(" Aligning...\n")    
    algs = []
    #prepare files - add timestamp
    tmpfn = os.path.join(tmpdir, "tmp.pid%s.%s"%(os.getpid(), time.time()))
    #write query
    tmpq = "%s.query" % tmpfn
    tmpqf = open(tmpq, "w")# as tmpqf:
    for i, (qid, qseq) in enumerate(zip(qids, qseqs), 1):
        tmpqf.write(">%s.%s\n%s\n"%(qid, i, qseq))
    tmpqf.close()
    #write targets
    tmpt = "%s.target" % tmpfn
    tmptf = open(tmpt, "w"); tmptf.write(matches); tmptf.close()
    #run blat
    tmpr = "%s.out" % tmpfn
    if dna:
        cmd  = "%s -q=dna -t=dna -noHead %s %s %s" % (blatpath, tmpt, tmpq, tmpr)
    else:
        cmd  = "%s -prot -noHead %s %s %s" % (blatpath, tmpt, tmpq, tmpr)
    if verbose:
        sys.stderr.write("  %s\n"%cmd)
    blatout = commands.getoutput(cmd)
    if verbose:
        sys.stderr.write("  %s\n"%blatout)
    #report error if not blat output
    if not os.path.isfile(tmpr):
        sys.stderr.write("BLAT failed (%s):\n%s\n"%(cmd, blatout))
        return
    #load results
    pq, algs = "", []
    for l in open(tmpr):
        ##BLAT PSL without header
        (matches, mismatches, repm, Ns, Qgapc, Qgaps, Tgapc, Tgaps, strand, \
         q, qsize, qstart, qend, t, tsize, tstart, tend, blocks, bsizes, \
         qstarts, tstarts) = l.split('\t')
        #unpack batch
        q = ".".join(q.split('.')[:-1])
        matches, mismatches = int(matches), int(mismatches)
        Tgapc, Tgaps = int(Tgapc), int(Tgaps)
        qstart, qend = int(qstart), int(qend)
        qstart, qend, qsize = int(qstart), int(qend), int(qsize)
        tstart, tend, tsize = int(tstart), int(tend), int(tsize)
        #get score, identity & overlap
        #score    = matches * 2 + mismatches * -1.5 + Tgapc * -11 + Tgaps * -1
        score    = matches * 5 + mismatches * -3 + Tgapc * -4 + Tgaps * -1
        alglen   = int(tend) - int(tstart)
        identity = 100.0 * matches / alglen
        overlap  = 100.0 * alglen / tsize
        #bitscore & evalue
        bitscore = (Lambda*score-math.log(K))/math.log(2)
        pvalue = evalue = 0
        if dblength:
            pvalue = 2**-bitscore
            evalue = len(qseqs[0]) * dblength * pvalue
        #get t and q ranges
        qranges  = get_ranges(qstarts, bsizes)
        tranges  = get_ranges(tstarts, bsizes)
        #if link provided, convert t into html link 
        if link:
            tlink = link % tuple([t]*link.count('%s'))
        else:
            tlink = t
        #yield if next query and any algs
        if q!=pq and algs:
            yield sorted(algs, key=lambda x: x[4], reverse=True)
            algs = []
        #add alg to list
        pq = q
        algs.append((q, tlink, round(identity,1), round(overlap,1), round(bitscore,2),\
                     "%.3g"%evalue, mismatches, Tgaps, alglen, qranges, tranges))
    #yield last alg
    if algs:
        yield sorted(algs, key=lambda x: x[4], reverse=True)
    #clean-up
    if "Error:" not in blatout and "No such file" not in blatout:
        os.unlink(tmpq)
        os.unlink(tmpt)
        os.unlink(tmpr)
        
def get_ranges(starts, sizes, offset=1):
    """Return str representation of alg ranges"""
    ranges = []
    for start, size in zip(starts.split(',')[:-1], sizes.split(',')[:-1]):
        start, size = int(start), int(size)
        start += offset
        end = start + size - offset
        coords = "%s-%s"%(start, end)
        ranges.append(coords)
    return " ".join(ranges)
    
def algs2formatted_output(algs, ifrac, ofrac, html, verbose, no_query):
    """Return TXT or HTML formatted output for matched sequences"""
    #get output
    blastTable = htmlTable()
    blastTable.style = "gtable"
    #add header
    headerNames = ["Query", "Target", "% identity", "% overlap", "Bit score", "E-value",\
                   "Mis- matches", "Gaps", "Alg. length", "Q. ranges", "T. ranges"]
    for name in headerNames: 
        blastTable.add_cell(0, name)
    #add results
    bident, boverlp = algs[0][2], algs[0][3]
    row_i = 0
    for hit in algs:
        if hit[2]<ifrac*bident or hit[3]<ofrac*boverlp:
            continue
        row_i += 1
        for item in hit:
            blastTable.add_cell(row_i, item)
    #rm query column
    if no_query:
        blastTable.remove_column(0)
    #return html or txt
    if html:
        return blastTable.asHTML()
    else:
        return blastTable.asTXT()

def fasta2hits(cur, db, table, seqcmd, qids, qseqs, blatpath, tmpdir, kmer, step, \
               seqlimit, html, link, sampling, dblength, verbose, dna=False, \
               ifrac=0.3, ofrac=0.3, dtype="uint32", no_query=True):
    """Report hits to qseq from hash table and sequences.
    Deals with both, single seq and list of translated sequences.
    """
    t0 = time.time()
    #deal with single sequence
    if type(qseqs) is str:
        qseqs = (qseqs, )
    #DNA or amino query
    if dna:
        seq2mers = dnaseq2mers
    else:
        seq2mers = aaseq2mers
    #get kmer matching sequences
    matches = seq2matches(cur, db, table, seqcmd, qids, qseqs, kmer, step, seqlimit, \
                          sampling, seq2mers, dtype, verbose)

    if not matches:
        return "#Your query (%s) didn't produce any hit.\n"%", ".join(qids)
    
    #align with blat
    out = ""
    algs = []
    for algs in hits2algs(qids, qseqs, matches, blatpath, tmpdir, link, dblength, dna, verbose):
        #return formatted output
        out += algs2formatted_output(algs, ifrac, ofrac, html, verbose, no_query)
    #write run stats
    dt = time.time() - t0
    info = "%s\t%s\t%s\t%.3f\t%s\n"%(datetime.ctime(datetime.now()), len("".join(qseqs)), len(algs), dt, "; ".join(qids))
    open(os.path.join(tmpdir, "fasta2hits.times.txt"), "a").write(info)
    if not out:
        return "#Your query (%s) didn't produce any valid hit.\n"%", ".join(qids)
    return out

def get_sixframe_translations(r):
    """Return list of six-frames translations"""
    seqs = []
    for seq in (r.seq, r.seq.reverse_complement()):
        for i in range(3): 
            ilen = 3 * ((len(seq)-i) // 3)
            seqs.append(str(seq[i:i+ilen].translate())) 
    return seqs

def get_random_sequence(cur, db, n, nprotids, seqcmd, verbose):
    """Return n random sequences"""
    if verbose:
        sys.stderr.write("Fetching %s random sequences...\n"%n)
    #get random ids
    if os.path.isfile(db):
        protids = random.sample(xrange(1, nprotids+1), n)
        seqs = sqlite2seq(cur, db, protids)
    else:
        #get no. of proteins in MySQL
        seqcmd = seqcmd.split(' where ')[0]
        s, e = seqcmd.upper().index('SELECT ')+7, seqcmd.upper().index(' FROM')
        cmd = seqcmd[:s] + "COUNT(*)" + seqcmd[e:]#; print cmd
        cur.execute(cmd)
        nprotids, = cur.fetchone()
        if verbose:
            sys.stderr.write(" %s proteins found...\n"%nprotids)
        #get random ids
        protids = random.sample(xrange(1, nprotids+1), n)
        seqs = []
        for x in protids:
            #print seqcmd+" LIMIT %s, 1"%x
            cur.execute(seqcmd+" LIMIT %s, 1"%x)
            seqs.append(">%s\n%s\n"%cur.fetchone())
        seqs = "".join(seqs)
    return seqs
        
def main():
    #compatible with Python2.6 without argparse
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)

    seqformats = ['fasta', 'fastq', 'gb', 'genbank', 'embl']
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='%(prog)s 1.4a')   
    parser.add_argument("-i", "--input",     type=file, default=sys.stdin, 
                        help="fasta stream    [stdin]")
    parser.add_argument("-b", "--batch",     type=int, default=1, 
                        help="query batch     [%(default)s]")
    parser.add_argument("--seqformat",       default="auto", choices = seqformats,
                        help="input format    [%(default)s]")
    parser.add_argument("-X", "--rapsiX",    default=False, action="store_true", 
                        help="6-frames translation of nucleotide query")
    parser.add_argument("--dna",             default=False, action='store_true',
                        help="DNA alphabet    [amino acids]")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("--html",            default=False, action="store_true", 
                        help="return HTML     [txt]")
    parser.add_argument("--no_query",        default=False, action="store_true", 
                        help="no query in output")
    parser.add_argument("--link",            default="", 
                        help="add html link matches [%(default)s]")
    parser.add_argument("-l", "--limit",     default=0, type=int,
                        help="stop after      [all]")
    parser.add_argument("--random",          default=0, type=int,
                        help="return N random sequence(s) and exit [%(default)s]")
    similo = parser.add_argument_group('Similarity search options')
    similo.add_argument("-k", "--kmer",      default=5, type=int, 
                        help="hash length     [%(default)s]")
    similo.add_argument("--blatpath",        default='blat', 
                        help="BLAT path       [%(default)s]")
    similo.add_argument("--tmpdir",          default='.', 
                        help="TEMP path       [%(default)s]")
    similo.add_argument("-s", "--step",       default=5, type=int, 
                        help="hash steps      [%(default)s]")
    similo.add_argument("--seqlimit",         default=100, type=int, 
                        help="max. seqs to retrieve [%(default)s]")
    similo.add_argument("-n", "--sampling",   nargs="*", default=[100], type=int,
                        help="sample n kmers per query %(default)s")
    similo.add_argument("--ifrac",            default=0.3, type=float, 
                        help="min. identity of best hit [%(default)s]")
    similo.add_argument("--ofrac",            default=0.3, type=float, 
                        help="min. overlap of best hit [%(default)s]")
    sqlopt = parser.add_argument_group('MySQL/SQLite options')
    sqlopt.add_argument("-d", "--db",         default="metaphors_201601", 
                        help="database        [%(default)s]")
    sqlopt.add_argument("-t", "--table",      default='hash2protids',
                        help="hashtable name  [%(default)s]")
    mysqlo = parser.add_argument_group('MySQL options')
    mysqlo.add_argument("-u", "--user",      default="script", 
                        help="database user   [%(default)s]")
    mysqlo.add_argument("-p", "--pswd",      default="python", 
                        help="user password   [%(default)s]")
    mysqlo.add_argument("-P", "--port",       default=4042, type=int, 
                        help="server port     [%(default)s]")
    mysqlo.add_argument("--host",            default='cgenomics.crg.es', 
                        help="database host   [%(default)s]")
    parser.add_argument("--seqcmd",           default="select protid, seq from protid2seq where protid in (%s)", 
                        help="SQL to fetch sequences [%(default)s]")
    mysqlo.add_argument("--dtype",            default="uint32", 
                        help="numpy data type for protids [%(default)s] ie. S8 for VARCHAR(8) or uint16 for SMALLINT UNSIGNED")
    mysqlo.add_argument("-z", "--dblength",   default=0, type=int,
                        help="database length for E-value calculation [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    if os.path.isfile(o.db):
        cnx = sqlite3.connect(o.db)
        #enable utf8 handling
        cnx.text_factory = str 
        cur = cnx.cursor()
        #get dblength & no. of proteins
        cur.execute("select value from meta_data where key='dblength'")
        o.dblength = int(cur.fetchone()[0])
        cur.execute("select value from meta_data where key='count'")
        proteins = int(cur.fetchone()[0])        
        if o.verbose:
            sys.stderr.write(" %s letters in %s entries\n"%(o.dblength, proteins))
    else:
        if not o.seqcmd:
            sys.stderr.write("Sequence command (--seqcmd) needed for MySQL\n")
            sys.exit(1)
        cnx = MySQLdb.connect(db=o.db, host=o.host, port=o.port, user=o.user, passwd=o.pswd)
        cur = cnx.cursor()
        proteins = 0

    #return random sequence(s)
    if o.random:
        o.output.write(get_random_sequence(cur, o.db, o.random, proteins, \
                                           o.seqcmd, o.verbose))
        return
        
    #handle gz/bz2
    handle = o.input
    seqformatpos = -1
    if handle.name.endswith('.gz'):
        seqformatpos = -2
        handle = gzip.open(handle.name)
    elif handle.name.endswith('.bz2'):
        import bz2
        seqformatpos = -2
        handle = bz2.BZ2File(handle.name)
    
    #get seqformat
    seqformat = o.seqformat
    if seqformat == "auto":
        #assume fasta on stdin
        if   handle == sys.stdin or handle.name.split('.')[seqformatpos]=='fa':
            seqformat = "fasta"
        elif handle.name.split('.')[seqformatpos]=='fq':
            seqformat = "fastq"
        elif handle.name.split('.')[seqformatpos] in seqformats:
            seqformat = handle.name.split('.')[seqformatpos]
        else:
            sys.exit("Cannot guess input sequence format: %s\n"%handle.name)
    
    #process entries        
    seqs, qids = [], []
    for i, r in enumerate(SeqIO.parse(handle, seqformat), 1):
        if not i%o.batch:
            sys.stderr.write(" %s \r"%i)
        if o.limit and i>o.limit:
            break
        #get 6-frames translations
        if o.rapsiX:
            seqs += get_sixframe_translations(r)
            qids += [r.id]*6
        else:
            seqs.append(str(r.seq))
            qids.append(r.id)
        #batch query
        if not i%o.batch:
            out = fasta2hits(cur, o.db, o.table, o.seqcmd, qids, seqs, o.blatpath, \
                             o.tmpdir, o.kmer, o.step, o.seqlimit, o.html, o.link, \
                             o.sampling, o.dblength, o.verbose, o.dna, \
                             o.ifrac, o.ofrac, o.dtype, o.no_query)
            o.output.write(out)
            seqs, qids = [], []
    if seqs:
        out = fasta2hits(cur, o.db, o.table, o.seqcmd, qids, seqs, o.blatpath, \
                         o.tmpdir, o.kmer, o.step, o.seqlimit, o.html, o.link, \
                         o.sampling, o.dblength, o.verbose, o.dna, \
                         o.ifrac, o.ofrac, o.dtype, o.no_query)
        o.output.write(out)
        
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("%s\n"%e)#\nI/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
