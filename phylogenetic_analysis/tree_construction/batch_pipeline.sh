#!/bin/sh

minimap2 -a -x asm20 --sam-hit-only --secondary=no --score-N=0 -t 10 \
        NC_045512.fas \
        ${name}.fasta \
        -o ${name}.sam

gofasta-linux-amd64 sam toMultiAlign \
       -s ${name}.sam \
       -t 10 \
       --reference NC_045512.fas \
       --trimstart 265 \
       --trimend 29674 \
       --trim \
       --pad > ${name}.aligned.fasta

python script/replace_N2hyphen.py \
       ${name}.aligned.fasta.filtered \
       > ${name}.aligned.N2hyphen.fasta

trimal \
  -in ${name}.aligned.N2hyphen.fasta \
  -out ${name}.aligned.N2hyphen.trimmed.fasta \
  -htmlout /dev/null \
  -gt 0.1

raxmlHPC-PTHREADS-SSE3 \
  -f a -x 12345 -p 12345 -# 3 -m GTRCAT -T 18 \
    -s ${name}.aligned.N2hyphen.trimmed.fasta \
    -n ${name}_tree

R --vanilla --slave --args \
  EPI_ISL_466615 \
  ${name}.aligned.N2hyphen.trimmed.fasta \
  RAxML_bipartitions.${name}_tree \
  ${name}.aligned.N2hyphen.trimmed.2nd.fasta \
  < ~/cov-phylo/omicron_tree/BA2/script/detect_outlier_OTU.R

raxmlHPC-PTHREADS-SSE3 \
  -f a -x 12345 -p 12345 -# 100 -m GTRCAT -T 18 \
    -s {name}.aligned.N2hyphen.trimmed.2nd.fasta \
    -n {name}.aligned.N2hyphen.trimmed.2nd_tree


