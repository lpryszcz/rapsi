#!/usr/bin/env python
desc="""Scan hash table and report matches. For searches for relatively similar
sequences, sequence sampling is recommended.

PARAMETERS EXAMPLES:
 --link  "<a target=_blank href='/?q=SeqInfo&seqid=Phy%s'>Phy%s</a>"

python wsgi/fasta2hits.py -vi ../test.fa --host cgenomics.crg.es -uscript -ppython --seqcmd "select concat(p.protid,'_',code), seq from protid2taxid p join protid2seq ps join species s on p.protid=ps.protid and p.taxid=s.taxid where p.protid in (%s)"

"""
epilog="""Author:
l.p.pryszcz@gmail.com

Barcelona/Mizerow, 13/11/2013

TO ADD:
- url input scanning for str threats like ", '


CHANGELOG:
v1.1:
- implemented nucleotide query (rapsiX)
- FastQ, genbank, embl support
- auto query format recognition

"""

import commands, gzip, os, sys, time, zlib
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
    
def seq2matches(cur, db, table, seqcmd, qid, qseqs, kmer, step, seqlimit, sampling, \
                intprotids, verbose):
    """Return matching protids and sequences"""
    if verbose:
        sys.stderr.write("Parsing %s aminos from: %s...\n" % (len("".join(qseqs)), qid))
    mers = set()
    for qseq in qseqs:
        mers = mers.union(seq2mers(qseq, kmer, step))
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
        if not seqcmd:
            merprotids = zlib.decompress(merprotids)
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
    
def hits2algs(qname, qseqs, matches, blatpath, tmpdir, link, verbose):
    """Align query with hits and return global algs."""
    if verbose:
        sys.stderr.write(" Aligning...\n")    
    algs = []
    #prepare files - add timestamp
    tmpfn = os.path.join(tmpdir, "tmp.%s"%time.time())
    #write query
    tmpq = "%s.query" % tmpfn
    with open(tmpq, "w") as tmpqf:
        for i, qseq in enumerate(qseqs, 1):
            tmpqf.write(">%s.%s\n%s\n"%(qname, i, qseq))
    #write targets
    tmpt = "%s.target" % tmpfn
    with open(tmpt, "w") as tmptf:
        tmptf.write("".join(matches))
    #run blat
    tmpr = "%s.out" % tmpfn
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

def fasta2hits(cur, db, table, seqcmd, qid, qseqs, blatpath, tmpdir, kmer, step, \
               seqlimit, html, link, sampling, intprotids, verbose):
    """Report hits to qseq from hash table and sequences.
    Deals with both, single seq and list of translated sequences.
    """
    t0 = time.time()
    #deal with single sequence
    if type(qseqs) is str:
        qseqs = (qseqs, )
    #get kmer matching sequences
    matches = seq2matches(cur, db, table, seqcmd, qid, qseqs, kmer, step, seqlimit, \
                          sampling, intprotids, verbose)

    if not matches:
        return "#Your query %s didn't produce any hit.\n"%qid
    
    #align with blat
    algs = hits2algs(qid, qseqs, matches, blatpath, tmpdir, link, verbose)

    if not algs:
        return "#Your query %s didn't produce any valid hit.\n"%qid
    
    #return formatted output
    out  = algs2formatted_output(algs, html, verbose)
    #write run stats
    dt = time.time() - t0
    with open(os.path.join(tmpdir, "fasta2hits.times.txt"), "a") as outtmp:
        outtmp.write("%s\t%s\t%s\n"%(qid, len(qseqs), dt))
    return out

def get_sixframe_translations(r):
    """Return list of six-frames translations"""
    seqs = []
    for seq in (r.seq, r.seq.reverse_complement()):
        for i in range(3):
            ilen = 3 * ((len(seq)-i) // 3)
            seqs.append(str(seq[i:i+ilen].translate())) 
    return seqs
    
def main():
    #compatible with severs without argparse
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)

    seqformats = ['fasta', 'fastq', 'gb', 'genbank', 'embl']
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='%(prog)s 1.1')   
    parser.add_argument("-i", "--input",     type=file, default=sys.stdin, 
                        help="fasta stream    [stdin]")
    parser.add_argument("--seqformat",       default="auto", choices = seqformats,
                        help="input format    [%(default)s]")
    parser.add_argument("-X", "--rapsiX",    default=False, action="store_true", 
                        help="6-frames translation of nucleotide query")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("--html",            default=False, action="store_true", 
                        help="return HTML     [txt]")
    parser.add_argument("--link",            default="", 
                        help="add html link matches [%(default)s]")
    parser.add_argument("-l", "--limit",     default=0, type=int,
                        help="stop after      [all]")
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
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    if os.path.isfile(o.db):
        cnx = sqlite3.connect(o.db)
        cur = cnx.cursor()
    else:
        if not o.seqcmd:
            sys.exit("Sequence command (--seqcmd) needed for MySQL")
        cnx = MySQLdb.connect(db=o.db, host=o.host, user=o.user, passwd=o.pswd)
        cur = cnx.cursor()

    #handle gz
    handle = o.input
    seqformatpos = -1
    if handle.name.endswith('.gz'):
        seqformatpos = -2
        handle = gzip.open(handle.name)
        
    #get seqformat
    seqformat = o.seqformat
    if seqformat == "auto":
        #assume fasta on stdin
        if   handle == sys.stdin or handle.name.split('.')[seqformatpos]=='.fa':
            seqformat = "fasta"
        elif handle.name.split('.')[seqformatpos]=='.fq':
            seqformat = "fastq"
        elif handle.name.split('.')[seqformatpos] in seqformats:
            seqformat = handle.name.split('.')[seqformatpos]
        else:
            sys.exit("Cannot guess input sequence format: %s\n"%handle.name)
    
    #process entries        
    for i, r in enumerate(SeqIO.parse(handle, seqformat), 1):
        if o.limit and i>o.limit:
            break
        #get 6-frames translations
        if o.rapsiX:
            seqs = get_sixframe_translations(r)
        else:
            seqs = str(r.seq)
        out = fasta2hits(cur, o.db, o.table, o.seqcmd, r.id, seqs, o.blatpath, \
                         o.tmpdir, o.kmer, o.step, o.seqlimit, o.html, o.link, \
                         o.sampling, o.intprotids, o.verbose)
        o.output.write(out)
        
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("\nI/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
