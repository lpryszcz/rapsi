#!/usr/bin/env python

#USAGE: cat pfam.fa.out.tsv | python pfam.benchmark.py Pfam-A.regions.tsv.gz

import gzip, os, sys

#C9WBN1	1	D19ACD554C340860	00d2a74435af6264b33f300212bbf72f	PF00001	57	331
#D8V2M4	1	695742C9F78E26B7	05a27745f5189888dbb43c4dc9477740	PF00001	47	105

def load_families(fn):
  """Return dict of gene2fams"""
  gene2fams = {}
  for i, l in enumerate(gzip.open(fn), 1):
    if not i%10000:
      sys.stderr.write(' %s  \r'%i)
    #if i>1000000: break
    ldata = l.split('\t')
    gene  = ldata[0]
    family = ldata[4]
    if gene in gene2fams:
      gene2fams[gene].add(family)
    else:
      gene2fams[gene] = set([family])
  print "%s gene2fams loaded." % len(gene2fams)
  return gene2fams

def get_match(matches, query):
  """Return first match different from query"""
  for match in matches:
    if match!=query:
      return match
  return None

def id2short(id):
  return id.split('.')[0]

def check_match(scores, query, hits):
  """Return updated scores"""
  #print query, id2short(query)
  m1 = get_match(hits[0], query)
  m2 = get_match(hits[1], query)
  if id2short(query) not in gene2fams or id2short(m1) not in gene2fams or id2short(m2) not in gene2fams: 
    return scores
  scores[0] += 1

  if gene2fams[id2short(query)].intersection(gene2fams[id2short(m1)]):
    scores[1] += 1

  if gene2fams[id2short(query)].intersection(gene2fams[id2short(m2)]):
    scores[2] += 1
  if m1 != m2:
    scores[3] += 1
  print scores    
  return scores

fn = sys.argv[1]

#load families
gene2fams = load_families(fn)

scores = [0, 0, 0, 0]
query = ""
hits = [[],[]]
for l in sys.stdin:
  ldata = l.split('\t')
  if l.startswith('#'):
    #store info
    if query:
      scores = check_match(scores, query, hits)
    hits = [[],[]]
    query = ""
  #rapsi
  elif len(ldata)==9:
    #Q4PA37.1        71.7    43.0    46      13      0       46      46-91   36-81
    hits[0].append(ldata[0])

  #blast  
  else:
    #F8Q7P8.1        F8Q7P8.1        100.00  114     0       0       1       114     1       114     1e-79    236
    hits[1].append(ldata[1])
    query = ldata[0]

if query:
  scores = check_match(scores, query, hits)
print scores
