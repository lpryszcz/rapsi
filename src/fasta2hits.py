#!/usr/bin/env python
desc="""Scan hash table and report matches.

For searches for relatively similar sequences, sequence sampling is recommended.

PARAMETERS EXAMPLES:
 --link  "<a target=_blank href='/?q=SeqInfo&seqid=Phy%s'>Phy%s</a>"

TO ADD:
- url input scanning for str threats like ", '

python wsgi/fasta2hits.py -vi ../test.fa --host cgenomics.crg.es -uscript -ppython --seqcmd "select concat(p.protid,'_',code), seq from protid2taxid p join protid2seq ps join species s on p.protid=ps.protid and p.taxid=s.taxid where p.protid in (%s)"
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona, 13/11/2013
"""

import commands, os, sys, time, zlib
import MySQLdb, resource, sqlite3
from datetime import datetime
from Bio import SeqIO
from htmlTable import htmlTable
from fasta2hash import seq2mers, b62decode
        
def sqlite2seq(cur, db, protids):
    """Return target fastas for protids from sqlite3."""
    #open target files
    cur.execute("SELECT name FROM file_data")
    files = {name: open(os.path.join(os.path.dirname(db), name)) for name, in cur.fetchall()}
    #get targets
    cmd = """SELECT f.name, offset, length FROM offset_data o JOIN file_data f
    ON o.file_number=f.file_number WHERE key IN (%s)""" % ",".join(str(p) for p in protids)
    cur.execute(cmd)
    targets = []
    for name, offset, length in cur.fetchall():
        files[name].seek(offset)
        targets.append(files[name].read(length))
    return "".join(targets)

def mysql2seq(cur, protids, seqcmd, intprotids):
    """Return target fasta for protids from mysql."""
    if intprotids:
        protidstext = ",".join(str(p) for p in protids)
    else:
        protidstext = "'" + "','".join(str(p) for p in protids) + "'"
    cur.execute(seqcmd % protidstext)
    return [">%s\n%s\n" % tup for tup in cur.fetchall()]
    
def seq2matches(cur, db, table, seqcmd, qid, qseq, kmer, step, seqlimit, sampling, \
                intprotids, verbose):
    """Return matching protids and sequences"""
    if verbose:
        sys.stderr.write("Parsing %s aminos from: %s...\n" % (len(qseq), qid))
    #mers = seq2mers(qseq, kmer, step, seqlimit, nosampling, verbose)
    mers = seq2mers(qseq, kmer, step)
    if verbose:
        sys.stderr.write("  %s mers: %s...\n"% (len(mers),", ".join(list(mers)[:6])))
    if sampling and len(mers)>sampling:
        mers = list(mers)[:sampling]
        if verbose:
            sys.stderr.write("  Sampled %s mers: %s...\n"% (len(mers),", ".join(list(mers)[:6])))
    #get protids for set of mers
    cmd = "select protids from `%s` where hash in ('%s')" % (table, "','".join(mers))    
    cur.execute(cmd)
    protids = {}
    for merprotids, in cur:
        #decode zlib coded str
        if not seqcmd: merprotids = zlib.decompress(merprotids)
        for protid in merprotids.split():
            #decode b62 int only for sqlite3
            #if not seqcmd: protid = b62decode(protid)
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

    if not fprotids:
        return 
    #get sequences
    if verbose:
        sys.stderr.write(" Fetching %s sequences...\n" % len(fprotids))
    #select protids as text or str
    if seqcmd:
        targets = mysql2seq(cur, fprotids, seqcmd, intprotids)
    else:
        targets = sqlite2seq(cur, db, fprotids)
    return targets
    
def hits2algs(qname, qseq, matches, blatpath, tmpdir, link, verbose):
    """Align query with hits and return global algs."""
    if verbose:
        sys.stderr.write(" Aligning...\n")    
    algs = []
    #prepare files - add timestamp
    tmpfn = os.path.join(tmpdir,"tmp.%s"%time.time())#pid%s"%os.getpid())
    tmpq = "%s.query" % tmpfn
    tmpqf = open(tmpq, "w")
    tmpqf.write(">%s\n%s\n"%(qname, qseq))
    tmpqf.close()
    tmpt = "%s.target" % tmpfn
    tmptf = open(tmpt, "w")
    tmptf.write("".join(matches))
    tmptf.close()
    tmpr = "%s.out" % tmpfn
    #run blat
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
    algs = []
    for l in open(tmpr):
        ##BLAT PSL
        (matches, mismatches, repm, Ns, Qgapc, Qgaps, Tgapc, Tgaps, strand, \
         q, qsize, qstart, qend, t, tsize, tstart, tend, blocks, bsizes, \
         qstarts, tstarts) = l.split('\t')
        matches, mismatches = int(matches), int(mismatches)
        Tgapc, Tgaps = int(Tgapc), int(Tgaps)
        qstart, qend = int(qstart), int(qend)
        qstart, qend, qsize = int(qstart), int(qend), int(qsize)
        tstart, tend, tsize = int(tstart), int(tend), int(tsize)
        #get score, identity & overlap - E-value???
        score    = matches * 2 + mismatches * -1.5 + Tgapc * -11 + Tgaps * -1
        e = '-'
        alglen   = int(tend) - int(tstart)
        identity = 100.0 * matches / alglen
        overlap  = 100.0 * alglen / tsize
        #get t and q ranges
        qranges  = get_ranges(qstarts, bsizes)
        tranges  = get_ranges(tstarts, bsizes)
        #if link provided, convert t into html link 
        if link:
            tlink = link % tuple([t]*link.count('%s'))
        else:
            tlink = t
        algs.append((tlink, round(identity, 1), round(overlap, 1), int(score), \
                     mismatches, Tgaps, alglen, qranges, tranges))
    #clean-up
    os.unlink(tmpq); os.unlink(tmpt); os.unlink(tmpr)
    return sorted(algs, key=lambda x: float(x[3]), reverse=True)

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
    
def algs2formatted_output(algs, html, verbose):
    """Return TXT or HTML formatted output for matched sequences"""
    #get output
    blastTable = htmlTable()
    blastTable.style = "gtable"
    #add header
    '''headerNames = ["Hit ID", "% identity", "% overlap", "E-value", "Score", \
                   "Alg. length", "Mis-matches", "Gaps", "Q. start", "Q. end", \
                   "T. start", "T. end"] '''
    headerNames = ["ID", "% identity", "% overlap", "Score", \
                   "Mis-matches", "Gaps", "Alg. length", "Q. ranges", "T. ranges"]
    for name in headerNames: 
        blastTable.add_cell(0, name)
    #add results
    bident, boverlp = algs[0][1], algs[0][2]
    row_i = 0
    for hit in algs:
        if hit[1]<0.3*bident or hit[2]<0.3*boverlp:
            continue
        row_i += 1
        for item in hit:
            blastTable.add_cell(row_i, item)
    #print
    #if verbose:
    #    sys.stderr.write(blastTable.asTXT())
    #return
    if html:
        return blastTable.asHTML()
    else:
        return blastTable.asTXT()

def fasta2hits(cur, db, table, seqcmd, qid, qseq, blatpath, tmpdir, kmer, step, \
               seqlimit, html, link, sampling, intprotids, verbose):
    """Report hits to qseq from hash table and sequences. """
    t0 = time.time()
    #get kmer matching sequences
    matches = seq2matches(cur, db, table, seqcmd, qid, qseq, kmer, step, seqlimit, \
                          sampling, intprotids, verbose)

    if not matches:
        return "#Your query %s didn't produce any hit.\n"%qid
    
    #align with blat
    algs = hits2algs(qid, qseq, matches, blatpath, tmpdir, link, verbose)

    if not algs:
        return "#Your query %s didn't produce any valid hit.\n"%qid
    
    #return formatted output
    out  = algs2formatted_output(algs, html, verbose)
    dt = time.time() - t0
    open(os.path.join(tmpdir, "fasta2hits.times.txt"), "a").write("%s\t%s\t%s\n"%(qid,len(qseq),dt))
    return out
        
def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", "--input",     type=file, default=sys.stdin, 
                        help="fasta file(s)   [stdin]")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("--html",            default=False, action="store_true", 
                        help="return HTML     [txt]")
    parser.add_argument("--link",            default="", 
                        help="add html link matches [%(default)s]")
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
    similo.add_argument("-n", "--sampling",   default=0, type=int,
                        help="sample n kmers from query [all]")
    sqlopt = parser.add_argument_group('MySQL/SQLite options')
    sqlopt.add_argument("-d", "--db",         default="metaphors_201310", 
                        help="database        [%(default)s]")
    sqlopt.add_argument("-t", "--table",      default='hash2protids4',
                        help="hashtable name  [%(default)s]")
    mysqlo = parser.add_argument_group('MySQL options')
    mysqlo.add_argument("-u", "--user",      default="lpryszcz", 
                        help="database user   [%(default)s]")
    mysqlo.add_argument("-p", "--pswd",      default="", 
                        help="user password   [%(default)s]")
    mysqlo.add_argument("--host",            default='localhost', 
                        help="database host   [%(default)s]")
    parser.add_argument("--seqcmd",           default="",
                        help="SQL to fetch sequences [%(default)s]")
    mysqlo.add_argument("--intprotids",       default=False,
                        help="protids are INT [STRING]") 
                         
    o = parser.parse_args()

    if os.path.isfile(o.db):
        cnx = sqlite3.connect(o.db)
        cur = cnx.cursor()
    else:
        if not o.seqcmd:
            sys.exit("Sequence command (--seqcmd) needed for MySQL")
        cnx = MySQLdb.connect(db=o.db, host=o.host, user=o.user, passwd=o.pswd)
        cur = cnx.cursor()

    for r in SeqIO.parse(o.input, 'fasta'):
        out = fasta2hits(cur, o.db, o.table, o.seqcmd, r.description, str(r.seq), o.blatpath, \
                         o.tmpdir, o.kmer, o.step, o.seqlimit, o.html, o.link, \
                         o.sampling, o.intprotids, o.verbose)
        o.output.write(out)
        
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
