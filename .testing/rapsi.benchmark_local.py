#!/usr/bin/env python

#USAGE: rapsi.benchmark_local.py pfam.fa

import os, sys, time
#sys.path.append('src')
import random, sqlite3

def time_query(cmd):
  t0 = time.time()
  os.system(cmd)
  dt = time.time() - t0
  return dt

def get_random_fasta(db):
  """Return fasta of random entry"""
  #connect to db
  cnx = sqlite3.connect(db+'.db3')
  cur = cnx.cursor()
  count, = cnx.execute("select value from meta_data where key='count'").fetchone()
  count = int(count)  
  cur.execute("SELECT name FROM file_data")
  files = {name: open(name) for name, in cur.fetchall()}
  seqcmd = "SELECT f.name, offset, length FROM offset_data o JOIN file_data f ON o.file_number=f.file_number WHERE key=%s"
  entryid = random.randint(1, count)
  name, offset, length = cnx.execute(seqcmd % entryid).fetchone()
  files[name].seek(offset)
  return files[name].read(length)  
  
query = "tmp.query.fa"
#db    = "phylomedb5/BlastDB.fasta"
db = sys.argv[1]
blastcmd = "blastall -pblastp -i %s -d %s -m8 -e1e-05 >> %s.out.tsv 2>> %s.stderr.txt" % (query, db, db, db)
rapsicmd = "src/fasta2hits.py -i %s -d %s.db3 --seqlimit 100 >> %s.out.tsv 2>> %s.stderr.txt" % (query, db, db, db)

out = open(db+'.times.tsv', 'a')
out.write("#protid\tseq len\trapsi\tblastp\n")
for i in xrange(1, 101): 
  #get random entry
  fasta  = get_random_fasta(db)
  protid = fasta[1:].split()[0]
  qseq = "".join(fasta.split()[1:]).replace('\n','')
  with open(query, 'w') as tmp:
    tmp.write(fasta)
  sys.stderr.write(" %s %s    \r" % (i, protid))
  
  dt1 = time_query(rapsicmd)
  dt2 = time_query(blastcmd)

  out.write("%s\t%s\t%s\t%s\n" % (protid, len(qseq), dt1, dt2))

