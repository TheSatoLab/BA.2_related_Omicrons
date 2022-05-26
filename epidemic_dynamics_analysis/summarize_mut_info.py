#!/usr/bin/env python

import sys, re
argvs = sys.argv

f = open(argvs[1])

f.__next__()

print("Id\tmut")
Id_mut_d = {}
for line in f:
  line = line.strip().split("\t")
  Id = line[2]
  mut_line = line[14].replace('(','').replace(')','')
  mut_l = mut_line.split(',')
  for mut in mut_l:
    print("\t".join([Id,mut]))
