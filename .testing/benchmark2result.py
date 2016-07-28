#!/usr/bin/env python
#Compute sensitivity and positivity from PFAM-A benchmark

import gzip, sys

#load families
prot2fam = {}
for i, l in enumerate(gzip.open("Pfam-A.regions.tsv.gz"), 1):
  if not i%10e3:
    sys.stderr.write(" %s \r"%i)
  protid, fam = l.split()[0], l.split()[4] 
  if protid not in prot2fam:
    prot2fam[protid] = set()
  prot2fam[protid].add(fam)
print "%s prot2fam loaded" % len(prot2fam)

def check_hits(hits, qname):
  """Count true/false positives"""
  #get best non-query hit
  i = 0
  while i < len(hits) and hits[i]==qname:
    i += 1
  if i >= len(hits):
    sys.stderr.write("No other match than query for: %s\n"%qname)
    return 0
  hit = hits[i]
  #get fam
  hfam = qfam = ""
  if hit.split('.')[0] in prot2fam:
    hfam = prot2fam[hit.split('.')[0]]
  if qfam.split('.')[0] in prot2fam:
    qfam = prot2fam[hit.split('.')[0]]
  if hfam.intersection(qfam):
    return 1
  return 0

#parse results
fn = "pfam.fa.out.tsv"
hits = [[], []]
qname = ""
i = k1 = k2 = 0
for l in open(fn):
  if l.startswith('#'):
    if qname:
      i += 1
      k1 += check_hits(hits[0], qname)
      k2 += check_hits(hits[1], qname)
      if not i%10:
        print i, k1, k2
    hits = [[], []]
    qname = ""
    continue
  ldata = l[:-1].split('\t')
  if len(ldata) == 9:
    hits[0].append(ldata[0])
  else:
    hits[1].append(ldata[1])
    qname = ldata[0]

if qname:
  i += 1
  k1 += check_hits(hits[0], qname)
  k2 += check_hits(hits[1], qname)

print i, k1, k2
'''
#ID     % identity      % overlap       Score   Mis-matches     Gaps    Alg. length     Q. ranges       T. ranges
F8P6P8.1        100.0   100.0   228     0       0       114     1-114   1-114
F8Q7P8.1        100.0   100.0   228     0       0       114     1-114   1-114
D6RKM8.1        79.2    90.6    135     22      0       106     9-114   7-112
B0CSP0.1        78.5    95.5    133     23      0       107     8-114   3-109
D8PV70.1        75.9    97.3    125     26      0       108     7-114   2-109
E2LV51.1        84.1    95.3    118     13      0       82      8-89    5-86
Q55R60.1        71.2    92.0    95      24      6       104     9-91 99-113     5-87 94-108
G4TGY8.1        72.7    86.1    93      26      1       99      8-21 24-107     4-17 19-102
E6R7V9.1        68.3    91.2    84      27      6       104     9-91 99-113     5-87 94-108
G7DT09.1        71.1    76.9    82      24      0       83      10-92   2-84
B8P9F5.1        81.7    59.4    81      11      0       60      55-114  41-100
E7A1S1.1        66.3    74.8    64      19      10      86      9-38 46-91      4-33 44-89
Q4PA37.1        71.7    43.0    46      13      0       46      46-91   36-81
F8Q7P8.1        F8Q7P8.1        100.00  114     0       0       1       114     1       114     1e-79    236
F8Q7P8.1        F8P6P8.1        100.00  114     0       0       1       114     1       114     1e-79    236
F8Q7P8.1        B0CSP0.1        78.50   107     23      0       8       114     3       109     2e-60    187
F8Q7P8.1        D6RKM8.1        76.58   111     26      0       4       114     2       112     5e-60    186
F8Q7P8.1        D8PV70.1        75.93   108     26      0       7       114     2       109     7e-59    183
F8Q7P8.1        G4TGY8.1        68.87   106     32      1       6       111     2       106     2e-49    160
F8Q7P8.1        Q55R60.1        68.81   109     33      1       6       114     2       109     3e-49    159
F8Q7P8.1        E2LV51.1        84.15   82      13      0       8       89      5       86      2e-48    156
F8Q7P8.1        E6R7V9.1        66.06   109     36      1       6       114     2       109     3e-47    154
F8Q7P8.1        B8P9F5.1        70.30   101     29      1       14      114     1       100     2e-44    146
F8Q7P8.1        G7DT09.1        65.09   106     34      1       9       114     1       103     2e-43    144
F8Q7P8.1        E7A1S1.1        56.36   110     45      1       7       113     2       111     1e-36    127
'''

