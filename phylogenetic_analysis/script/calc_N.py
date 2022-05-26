#!/usr/bin/env python

import sys
argvs = sys.argv

f_name = argvs[1]
f = open(f_name)
start = int(argvs[2])
end = int(argvs[3])
N_limit = float(argvs[4])

out_fasta = open(f_name + ".filtered","w")
out_txt = open(f_name + ".txt","w")

len_interest = end - start

seq_d = {}
for line in f:
  line = line.strip()
  if '>' in line:
    name = line.replace('>','')
    seq_d[name] = ""
  else:
    seq_d[name] += line

for name in seq_d:
  seq = seq_d[name]
  count_N_interest = seq.count("N",start,end)
  prop_N_interest = float(count_N_interest) / len_interest
  res_l = [name,count_N_interest,prop_N_interest]
  print("\t".join([str(c) for c in res_l]),file=out_txt)
  if prop_N_interest < N_limit:
    print("\n".join([">" + name,seq]),file = out_fasta)
 
