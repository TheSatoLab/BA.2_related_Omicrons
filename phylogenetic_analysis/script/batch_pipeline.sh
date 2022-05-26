#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 18
#$ -l exclusive,s_vmem=10G,mem_req=10G
##$ -o /dev/null
##$ -e /dev/null
##$ -t 1-2:1


cd /home/jampei/cov-phylo/BA2_L452/2022_04_23/2nd/3rd/4th

#path
export PATH=${PATH}:~/bin
export PATH=${PATH}:/usr/local/package/fasttree/2.1.10/bin
export PATH=${PATH}:/home/satolab/trimAl/source
export PATH=${PATH}:/usr/local/package/minimap/2.17/bin
export PATH=${PATH}:/usr/local/package/raxml/8.2.12/bin
export OMP_NUM_THREADS=18


#name=BA2.with_BA1_BA3_B11.BA2_hap_sampled
#name=$(cat script/sample_list.txt | awk -v i=$SGE_TASK_ID 'BEGIN{FS="\t"}NR==i{print $1}')

<<COMMENTOUT
python script/fiter_tree_for_download_package.py \
  metadata.${name}.txt \
  ~/cov-phylo/gisaid/2022_04_23/sequences.fasta \
  > ${name}.fasta


#cat metadata.representative.BA1_BA2.fas \
#    metadata.chimera_candidates.fasta \
#    > metadata.chimera_candidates.catted.fasta


##COMMENTOUT
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


#mafft --thread 10 --6merpair --addfragments \
#      /home/jampei/cov-phylo/gisaid/2021-09-07_13-20.filtered/sequences.2021-09-07_13-20.filtered.delta.fasta \
#      NC_045512.fas \
#      > sequences.2021-09-07_13-20.filtered.delta_with_Wuhan.aligned.fasta


python script/calc_N.py \
       ${name}.aligned.fasta \
       21562 25384 0.02


python script/replace_N2hyphen.py \
       ${name}.aligned.fasta.filtered \
       > ${name}.aligned.N2hyphen.fasta
##COMMENTOUT

trimal \
  -in ${name}.aligned.N2hyphen.fasta \
  -out ${name}.aligned.N2hyphen.trimmed.fasta \
  -htmlout /dev/null \
  -gt 0.1

#COMMENTOUT
raxmlHPC-PTHREADS-SSE3 \
  -f a -x 12345 -p 12345 -# 3 -m GTRCAT -T 18 \
    -s ${name}.aligned.N2hyphen.trimmed.fasta \
    -n ${name}_tree

#COMMENTOUT

##R --vanilla --slave --args \
##  ${name}.aligned.N2hyphen.trimmed.fasta \
##  RAxML_bipartitions.${name}_tree \
##  ${name}.aligned.N2hyphen.trimmed.2nd.fasta \
##  < script/detect_outlier_OTU.R

##COMMENTOUT

R --vanilla --slave --args \
  EPI_ISL_466615 \
  ${name}.aligned.N2hyphen.trimmed.fasta \
  RAxML_bipartitions.${name}_tree \
  ${name}.aligned.N2hyphen.trimmed.2nd.fasta \
  < ~/cov-phylo/omicron_tree/BA2/script/detect_outlier_OTU.R

COMMENTOUT



##2nd tree construction######################
raxmlHPC-PTHREADS-SSE3 \
  -f a -x 12345 -p 12345 -# 3 -m GTRCAT -T 18 \
    -s BA2.with_BA1_BA3_B11.BA2_hap_sampled.aligned.N2hyphen.trimmed.2nd.rm_BA2_with_69_70del.fasta \
    -n BA2.with_BA1_BA3_B11.BA2_hap_sampled.aligned.N2hyphen.trimmed.2nd.rm_BA2_with_69_70del.tree







#FastTreeMP -nt -gtr \
#  clade_GRA.filtered.aligned.N2hyphen.trimmed.fasta \
#  > clade_GRA.filtered.aligned.N2hyphen.trimmed.nwk

