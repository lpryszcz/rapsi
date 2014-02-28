#!/usr/bin/env python

import httplib, sys, time, urllib
sys.path.append('../cgenomics/beta.phylomedb.org/wsgi/')
from rapsi import *

def http_query(fasta, conn, methodurl):
  """Return dt and server repsonse"""
  t0 = time.time()
  conn.request("GET", methodurl+urllib.quote(fasta))
  r = conn.getresponse()
  html = r.read() 
  if r.status != 202:
    time.sleep(5)
    return "error", r
  dt = time.time() - t0
  return dt, r

#connect
p = PhylomeDBConnector(host="cgenomics.crg.es", user="phyReader", passwd="phyd10.-Reader", db = "phylomedb_5")
cursor = p._SQLconnection.cursor()

#change params
url = "beta.phylomedb.org"
rapsiurl = "/wsgi/rapsi.py/similarity?fasta="
blasturl = "/wsgi/phylomedb.py/blast_hits?sequence="

conn = httplib.HTTPConnection(url)

out = open(sys.argv[0]+'.txt', 'a')
out.write("#protid\tseq len\trapsi\tblastp\n")
for i in xrange(1, 1001): 
  #get random entry
  metaid = p.get_random_seed(PUBLIC_PHYLOMES)
  seqinfo = p.get_seqid_info(metaid)
  protid, qseq = seqinfo['protid'], seqinfo['seq']
  fasta = ">%s\n%s" % (protid, qseq)
  sys.stderr.write(" %s %s    \r" % (i, protid))
  
  dt1, r1 = http_query(fasta, conn, rapsiurl)
  dt2, r2 = http_query(fasta, conn, blasturl)

  out.write("%s\t%s\t%s\t%s\n" % (protid, len(qseq), dt1, dt2))

