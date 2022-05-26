#!/usr/bin/env python

import sys, re
argvs = sys.argv

metadata_f = open(argvs[1])
fasta_f = open(argvs[2])

metadata_f.__next__()
name_Id_d = {}
for line in metadata_f:
  line = line.strip().split("\t")
  name = line[0]
  Id = line[2]
  name_Id_d[name] = Id

name_s = set(name_Id_d.keys())

seq_d = {}
frag = 0
for line in fasta_f:
  line = line.strip()
  if '>' in line:
    name = line.replace('>','')
    name = name.split('|')[0]
    frag = 0
    if name in name_s:
      frag = 1
  if frag == 1:
    if name not in seq_d:
      seq_d[name] = ""
      continue
    seq = line
    seq_d[name] += seq

for name in seq_d:
  if name in name_Id_d:
    Id = name_Id_d[name]
    seq = seq_d[name]
    res = ">" + Id + "\n" + seq
    print(res)
