#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(ggplot2)
library(rbin)
library(cmdstanr)
library(patchwork)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggtree)
library(ggnewscale)
library(ape)
library(castor)
library(phangorn)
library(igraph)
library(TreeTools)




args = commandArgs(trailingOnly=T)

mut.interest.v <- c("Spike_L452R","Spike_L452Q","Spike_L452M")

tree.name <- 'RAxML_bipartitions.BA2.with_outgroups.aligned.N2hyphen.trimmed.2nd.rerooted'
tree <- read.tree(tree.name)


tree.info.df <- ggtree(tree)$data

count_descendants <- function(node){
  n.descendants <- length(Descendants(tree, node, type = "tips")[[1]])
  return(n.descendants)
}

tree.info.df$num.decendants <- as.numeric(map(tree.info.df$node,count_descendants))

#k <- 10
extract_ancestral_nodes <- function(node){
  node.v <- Ancestors(tree, node, type = "all") #[1:k]
  node.v <- as.numeric(na.omit(node.v))
  return(node.v)
}

ancestral.nodes.l <- map(tree.info.df$node,extract_ancestral_nodes)


tree.tip.df <- data.frame(tip_Id = 1:length(tree$tip.label), Id = tree$tip.label)

n.seq <- length(tree$tip.label)

metadata.name <- 'metadata.tsv'   #paste('/media/jampei/C6700DA7700D9F75/data/Sato_analysis/BA2_L452/tree/2022_04_23/2nd/metadata.',name,'.2022_04_23.txt',sep="")
metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)
#metadata <- metadata %>% select(-mut)
Id.analyzed.v <- as.character(tree.tip.df$Id)


#mutation data
mut.info.name <- 'metadata.mut_long.tsv'
mut.info <- fread(mut.info.name,header=T,sep="\t",check.names=T)
mut.info <- mut.info %>% filter(Id %in% Id.analyzed.v)

mut.info.interest <- mut.info %>% filter(str_detect(mut,"Spike_L452"))

mut.info.interest.spread <- mut.info.interest %>% filter(mut %in% mut.interest.v) %>% mutate(value = 1) %>% spread(key = mut, value = value)

metadata <- metadata[match(tree.tip.df$Id,Accession.ID),]

metadata <- metadata %>% left_join(mut.info.interest.spread %>% rename(Accession.ID = Id),by="Accession.ID")
metadata[is.na(metadata)] <- 0

metadata <- metadata %>%
                       mutate(Collection.date = as.Date(Collection.date),
                              region = str_split(Location," / ",simplify = T)[,1],
                              country = str_split(Location," / ",simplify = T)[,2])



mut.name <- "Spike_L452R"

min.num.decendants <- 10

sum.df <- data.frame()

for(mut.name in mut.interest.v){

state.v <- metadata %>% pull(mut.name)

fit.Mk <- asr_mk_model(tree, state.v + 1, rate_model = "ARD")

state.mat <- fit.Mk$ancestral_likelihoods
state.node.v <- state.mat[,1]
state.node.v <- ifelse(state.node.v>0.5,0,1)

state.color <- c(state.v,state.node.v)
state.color <- factor(state.color,levels=c(0,1))




tree.info.df.interest <- tree.info.df %>% select(parent,node,num.decendants) %>% mutate(state = state.color)



tree.info.df.interest.parent <- tree.info.df.interest %>% select(node,state) %>% rename(parent = node, state.parent = state)
tree.info.df.interest.merged <- merge(tree.info.df.interest,tree.info.df.interest.parent,by="parent")
tree.info.df.interest.merged.0to1 <- tree.info.df.interest.merged %>% filter(state == 1, state.parent == 0, num.decendants >= min.num.decendants)


node.0to1.v <- NA
if(nrow(tree.info.df.interest.merged.0to1) >0){
  node.0to1.v <- tree.info.df.interest.merged.0to1$node
}

if(!is.na(node.0to1.v)){
for(j in 1:length(node.0to1.v)){
  node.0to1 <- node.0to1.v[j]
  ancestral.nodes.v <- ancestral.nodes.l[[node.0to1]][1:5]
  tree.info.df.interest.ancestors_1 <- tree.info.df.interest %>% filter(node %in% ancestral.nodes.v,state == 1)
  if(nrow(tree.info.df.interest.ancestors_1) > 0){
    node.0to1.v[j] <- NA
  }
}

node.0to1.v <- as.vector(na.omit(node.0to1.v))

if(length(node.0to1.v)>1){
  ancestral.nodes.k.l <- list()
  for(j in 1:length(node.0to1.v)){
    #j <- 1
    node.0to1 <- node.0to1.v[j]
    ancestral.nodes.v <- ancestral.nodes.l[[node.0to1]][1:5]
    ancestral.nodes.k.l[[j]] <- ancestral.nodes.v
  }

  link.df <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
  for(j1 in 1:length(ancestral.nodes.k.l)){
    for(j2 in 2:length(ancestral.nodes.k.l)){
      if(j1 < j2){
        overlap.v <- intersect(ancestral.nodes.k.l[[j1]],ancestral.nodes.k.l[[j2]])
        if(length(overlap.v)>0){
          temp.df <- data.frame(j1 = j1, j2 = j2)
          link.df <- rbind(link.df,temp.df)
        }

      }

    }
  }
  if(nrow(link.df)>0){
  g <- graph.data.frame(link.df,directed=F)
  dg.l <- decompose(g, min.vertices = 1)

  network_node.v.new <- c(NA)
  for(j in 1:length(dg.l)){
    dg <- dg.l[[j]]
    dg.info <- as.data.frame(degree(dg))
    network_node.v <- as.numeric(rownames(dg.info))
    node.0to1.v.clustered <- node.0to1.v[network_node.v]
    node.0to1.new <- getMRCA(tree,node.0to1.v.clustered)
    network_node.v.new <- c(network_node.v.new,node.0to1.new)
  }

  node.0to1.v <- network_node.v.new
  }
}



root_node <- 8112
root_node.state <- tree.info.df.interest[tree.info.df.interest$node == root_node,]$state %>% as.character() %>% as.numeric()
if(root_node.state == 1){
  node.0to1.v <- root_node
}

temp.df <- data.frame(mut.name = mut.name, node.0to1 = node.0to1.v)
temp.df <- na.omit(temp.df)

sum.df <- rbind(sum.df,temp.df)

state.color2 <- as.numeric(as.character(state.color))

if(nrow(temp.df) > 0){
  state.color2[temp.df$node.0to1] <- 2
}

state.color2 <- factor(state.color2,levels=c(0,1,2))


g <- ggtree(tree) + aes(color=state.color2)
g <- g + scale_color_manual(breaks=c(0,1,2),values = c("gray70","gold2","red"))
g

pdf.name <- paste(mut.name,'.pdf',sep="")
pdf(pdf.name,width=5,height=5)
plot(g)
dev.off()



}

}



sum.decendants.df <- data.frame()

for(i in 1:nrow(sum.df)){

mut.name <- sum.df[i,]$mut.name
node.interest <- sum.df[i,]$node.0to1

decendants.v <- Descendants(tree, node.interest, type = "tips")[[1]]
decendants.Id.v <- tree.tip.df %>% filter(tip_Id %in% decendants.v) %>% pull(Id)

temp.df <- data.frame(mut.name = mut.name, node.interest = node.interest, Id = decendants.Id.v)
sum.decendants.df <- rbind(sum.decendants.df,temp.df)
}


sum.decendants.df <- sum.decendants.df %>% inner_join(metadata %>% select(Id = Accession.ID, Pango.lineage, country, contains("Spike")),by="Id")





out.name <- 'sum.L452_gain_node.txt'
write.table(sum.decendants.df,out.name,col.names=T,row.names=F,sep="\t",quote=F)




