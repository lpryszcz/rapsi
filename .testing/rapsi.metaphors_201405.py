#!/usr/bin/env python

import commands, httplib, os, sys, time, urllib
#from rapsi import *

def http_query(q, conn, methodurl):
  """Return dt and server repsonse"""
  t0 = time.time()
  conn.request("GET", methodurl+urllib.quote(q))
  r = conn.getresponse()
  html = r.read() 
  if r.status != 202:
    time.sleep(5)
    return "error", r
  dt = time.time() - t0
  return dt, html

#change params
url = "betaorthology.phylomedb.org"
sequrl = "/wsgi/query.py/getFasta?metaid="

tmpfn = 'tmp_pid%s' % os.getpid()

fn = sys.argv[1]
ofn = fn+'.txt'
if os.path.isfile(ofn):
  sys.exit("File exists: %s!"%ofn)
out = open(ofn, 'w')
#overwrite blast/rapsi output
with open(tmpfn+'_blast_output', 'w') as tmp: tmp.write("#blast all vs metaphors\n")
with open(tmpfn+'_rapsi_output', 'w') as tmp: tmp.write("#rapsi vs metaphors\n")
for i, line in enumerate(open(fn), 1): 
  #get query
  query = line[:-1].split('\t')[4]
  if not query.startswith(">"):
    sys.stderr.write("[WARNING] Wrong query: %s\n"%query)
    continue
  sys.stderr.write(" %s %s    \n" % (i, query))
  metaid = query[1:].split('_')[0]
  conn = httplib.HTTPConnection(url)
  dt0, fasta = http_query(metaid, conn, sequrl)
  #write fasta
  with open(tmpfn+"_query", 'w') as tmp: tmp.write(fasta)
  #run rapsi
  cmd1="src/fasta2hit.py -d metaphors.fa.gz.db3 -i %s_query >> %s_rapsi_output"%(tmpfn, tmpfn)
  t0 = time.time()
  blastout = commands.getoutput(cmd1)
  dt1 = time.time() - t0  
  #run blast
  cmd2="blastall -p blastp -e 1e-05 -d metaphors -i %s_query -m 8 >> %s_blast_output"%(tmpfn, tmpfn)
  t0 = time.time()
  blastout = commands.getoutput(cmd2)
  dt2 = time.time() - t0
  out.write("%s\t%s\n" % (line[:-1], dt2))

