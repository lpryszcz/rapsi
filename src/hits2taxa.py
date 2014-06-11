#!/usr/bin/env python
desc="""Parse alignments (SAM or blast8) and report organism summary.
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 28/02/2014
"""

import argparse, gzip, os, re, sys
import numpy as np
import sqlite3
from datetime import datetime
from Bio      import SeqIO
from taxonomy import Taxonomy

#cigar patterns
spat = re.compile('\d+[HSI]') #split at hardclip/softclip/indels
cpat = re.compile('[MD]') #

def cigar2mlen(cigar):
    """Return match length from cigar. This is the last aligned
    ref position. May need adjustment!
    """
    #remove soft clipped and insertions into ref and the get
    aligned = cpat.split("".join(spat.split(cigar)))
    return sum(int(x) for x in aligned if x)

def get_matches_sam(samstream, verbose):
    """Parse SAM and yield matches"""
    pread = ""
    pmapq = 0
    hits  = []
    for l in samstream:
        #skip SAM header
        if l.startswith('@'):
            continue
        #rname, flag, ref, pos, mapq, cigar, mref, mpos, isize, seq, quals = l.split('\t')
        try:
            rname, flag, ref, pos, mapq, cigar = l.split('\t')[:6]
        except:
            sys.stderr.write("[Warning] Cannot parse SAM line: %s\n"%str(l.split('\t')[:6]))
            continue
        #you may also need comment (if spaces in read name)
        #unmapped
        if ref == "*":
            #use megablast optionally
            #just to count properly reads!
            yield None, []
            continue
        #
        flag, pos, mapq = int(flag), int(pos), int(mapq)
        #store
        if pread != rname:
            if hits:
                yield pread, hits 
                hits  = []
            pread = rname
            pmapq = 0
        #skip hits with lower mapping
        #consider some more sophisticated way of checking
        #as for long reads parts can align to many chr
        if mapq < pmapq:
            continue
        #get gi int
        if ref.startswith("gi|"):
            ref = int(ref.split("|")[1])
        #store hit and mapq
        pmapq = mapq
        mlen  = cigar2mlen(cigar)
        hits.append((ref, pos, pos+mlen))
            
    #store last bit
    if hits:
        yield pread, hits
    
def get_matches_blast8(input, verbose):
    """Return matches from blast8 or rapsi output."""
    pscore = q = pq = 0
    hits   = []    
    for line in input:
        if line.startswith("#") or q!=pq:
            #report previous hits
            if hits:
                yield pq, hits
                hits = []
                pscore = 0
            pq = q
            continue
        ldata = line[:-1].split('\t')
        #blast
        if len(ldata)==12:
            q, ref, ident, m, matches, gaps, qs, qe, ts, te, e, score = ldata
            start, end = int(ts), int(te)
            if ts>te:
                start, end = end, start
        #rapsi
        else:
            q, ref, ident, overlap, score, e, mmatches, gaps, alen, qranges, tranges = ldata
            start = int(tranges.split()[0].split('-')[0])
            end   = int(tranges.split()[-1].split('-')[1])
        if float(score) < pscore:
            continue
        #get gi int
        if ref.startswith("gi|"):
            ref = int(ref.split("|")[1])
        #store match
        hits.append((ref, start, end))
        pscore = float(score)
    if hits:
        yield "None", hits
     
def get_taxa(hits, taxa, verbose):
    """Return taxa and genes matched by given read.
    Return None if htaxid appears in hits.
    """
    taxid2gi = {}
    #process each match
    for gi, pos, end in hits:
        #get taxid
        taxid = taxa.gi2taxid(gi)
        if not taxid: 
            if verbose:
                sys.stderr.write("[Warning] No taxid for gi: %s\n"%gi)
            continue
        #add taxa
        mdata = (gi, pos, end)
        if taxid not in taxid2gi:
            taxid2gi[taxid] = [mdata]
        else:
            taxid2gi[taxid].append(mdata)
            
    #get most occuring taxa as representative
    #taxid = sorted(taxid2gi.keys(), key=lambda x: len(taxid2gi[x]), reverse=True)[0]

    #no taxa matched
    if   not taxid2gi:
        return None, []
    #one taxa matched - spp or strain
    elif len(taxid2gi) == 1:
        taxid, matches = taxid2gi.popitem()
    ##or genus or family instead of strain/species!
    else:
        taxid = taxa.collapse_taxa(taxid2gi.keys())
        if not taxid:
            return None, []
        matches = [mdata for t, mdata in taxid2gi.iteritems()]
    return taxid, matches
    
def hits2taxa(input, out, db, verbose, limit=0): 
    """Process fastq from input file.
    You may play with bufsize, so processes runs without waiting.
    """
    #init taxonomy
    taxa = Taxonomy(db)

    #hadle gzipped/bzip2 stream
    if input.name.endswith('.gz'):
        input = gzip.open(input.name)
    elif input.name.endswith('.bz2'):
        import bz2
        input = bz2.BZ2File(input.name)
    #get match generator
    if input==sys.stdin:
        line0 = input.readline()
        if line0.startswith('@'):
            mGenerator = get_matches_sam(input, verbose)
        else:
            mGenerator = get_matches_blast8(input, verbose)
    #get sam stream
    elif input.name.endswith(('.sam','.sam.gz')):
        mGenerator = get_matches_sam(input, verbose)
    else:
        mGenerator = get_matches_blast8(input, verbose)

    #process reads in 1K batches?print 
    if verbose:
        sys.stderr.write("[%s] Processing reads from %s ...\n"%(datetime.ctime(datetime.now()), input.name))
    #get taxa and genes
    taxid2reads   = {}
    taxid2matches = {}
    k = 0
    for i, (rname, hits) in enumerate(mGenerator, 1):
        if limit and i>limit:
            break
        if not rname:
            continue
        #print info
        if verbose and i%1e4 == 1:
            sys.stderr.write(" %s parsed. %.2f%s with taxa  \r"%(i, k*100.0/i, '%'))
        #get taxa
        taxid, matches = get_taxa(hits, taxa, verbose)
        if not taxid: 
            continue
        k += 1
        if taxid not in taxid2reads:
            taxid2reads[taxid] = 0 
        #store read name & genes
        taxid2reads[taxid] += 1

    #report
    if not taxid2reads:
        sys.exit("No matches found!")
    ##foreign reads
    freads = sum(reads for taxid, reads in taxid2reads.iteritems())
    header = "#name\ttaxid\treads\t%\n"
    out.write(header)
    out.write("%s\t%s\t%s\t%.2f\n"%("unknown", "-", i-freads, 100.0*(i-freads)/i))
    for taxid, reads in sorted(taxid2reads.iteritems(), key=lambda x: x[1], reverse=True)[:10]: 
        out.write("%s\t%s\t%s\t%.2f\n"%(taxa[taxid][1], taxid, reads, 100.0*reads/i))
    #print summary
    sys.stderr.write("[hits2taxa] %s entries processed!\n" % (i,))

def main():
    usage   = "src/%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)
  
    parser.add_argument("-v", dest="verbose",   default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-d", "--db",        default="taxonomy.db3",
                        help="taxonomy path  [%(default)s]")
    parser.add_argument("-i", "--input",     type=file, default=sys.stdin,
                        help="input sam stream  [stdin]")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType("w"), 
                        help="output stream     [stdout]")
    parser.add_argument("-l", "--limit",     default=0, type=int, 
                        help="no. of reads to process and rows of tables to load [all]")
                        
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("[hits2taxa] Options: %s\n"%str(o))

    hits2taxa(o.input, o.output, o.db, o.verbose, o.limit)
	
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\n[hits2taxa] Ctrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write( "[hits2taxa] Time elapsed: %s\n" % dt )
