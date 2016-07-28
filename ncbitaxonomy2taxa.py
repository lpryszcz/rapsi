#!/usr/bin/env python
desc="""Parse NCBI taxdump and report taxonomy table.

###
#Recommended usage (~20m):
sqlite3 taxonomy.db3 "create table taxa_info (taxid INTEGER PRIMARY KEY, parent_id INTEGER, name TEXT, rank TEXT)"
src/ncbitaxonomy2taxa.py -v | sqlite3 -separator $'\\t' taxonomy.db3 '.import /dev/stdin taxa_info' 
sqlite3 taxonomy.db3 "create table gi2taxid  (gi INTEGER PRIMARY KEY, taxid INTEGER)"
wget -O- ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_????.dmp.gz | zcat | sqlite3 -separator $'\\t' taxonomy.db3 '.import /dev/stdin gi2taxid'
###
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 28/02/2014
"""

import argparse, gzip, os, sys
from datetime import datetime
import tarfile
from string import strip

def ncbitaxonomy2taxa(taxdump, out, verbose):
    """Parse taxdump for NCBI and report taxa.gz."""
    #fetch taxdump.tar.gz
    if not os.path.isfile(os.path.basename(taxdump)):
        os.system('wget %s'%taxdump)
    #store info
    taxid2data = {}
    #open tarfile
    tar = tarfile.open(os.path.basename(taxdump))
    #get nodes
    if verbose:
        sys.stderr.write("Parsing nodes...\n")
    m = tar.getmember('nodes.dmp')
    for i, line in enumerate(tar.extractfile(m), 1):
        if verbose and not i%10000:
            sys.stderr.write(" %s \r"%i)
        taxid, ptaxid, rank =  map(strip, line.split("|"))[:3]
        taxid, ptaxid = int(taxid), int(ptaxid)
        taxid2data[taxid] = [ptaxid, rank]
    #get names
    m = tar.getmember('names.dmp')
    if verbose:
        sys.stderr.write("Parsing names...\n")
    for i, line in enumerate(tar.extractfile(m), 1):
        if verbose and not i%10000:
            sys.stderr.write(" %s \r"%i)
        fields =  map(strip, line.split("|"))        
        if fields[3] == "scientific name":
            taxid = int(fields[0])
            name  = fields[1]
            ptaxid, rank = taxid2data[taxid]
            out.write("%s\t%s\t%s\t%s\n"%(taxid, ptaxid, name, rank))
    #close output
    out.close()

def main():
    usage  = "src/%(prog)s -v"
    parser = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                     formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", dest="verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", "--input",  default="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", 
                        help="input stream [%(default)s]")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("w"), 
                        help="output stream [stdout]")
                         
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % str(o))

    ncbitaxonomy2taxa(o.input, o.output, o.verbose)
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\n Ctrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
            
    
