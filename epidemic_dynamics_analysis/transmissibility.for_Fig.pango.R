#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(ggplot2)
library(rbin)
library(cmdstanr)
library(patchwork)
library(RColorBrewer)

args = commandArgs(trailingOnly=T)

##########args##########
#input
download.date <- "2022-05-15"
download.date <- as.Date(download.date)
stan_f.name <- 'multinomial_independent.stan' #args[2]
metadata.name <- 'metadata.tsv' #args[4]
mut.info.name <- 'metadata.mut_long.tsv'
mut.interest <- "Spike_L452"


#output
out.prefix <- '2022_05_15_Fig.pango' #args[5]
pdf.observed.name <- paste(out.prefix,".method1.observed2.pdf",sep="")
pdf.theta.name <- paste(out.prefix,".method1.theta2.pdf",sep="")
pdf.growth.rate.name <- paste(out.prefix,".method1.growth_rate2.pdf",sep="")
txt.growth.rate.name <- paste(out.prefix,".method1.growth_rate2.txt",sep="")


##########parameters##########
#general
core.num <- 4
variant.ref <- "BA.2"

#period to be analyzed
day.delay <- 0
day.analyzed <- 100
day.recent <- 14

date.end <- download.date - day.delay
date.start <- date.end - day.analyzed + 1
date.recent.start <- date.end - day.recent + 1


metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)
metadata.filtered <- metadata %>%
                       distinct(Accession.ID,.keep_all=T) %>%
                       filter(Host == "Human",
                              str_length(Collection.date) == 10,
                              #Is.low.coverage.!="True",
                              Pango.lineage != "",
                              Pango.lineage != "None",
                              Sequence.length > 29000,
                              !str_detect(Additional.location.information,"[Qq]uarantine")
                              )
 
metadata.filtered <- metadata.filtered %>%
                       mutate(Collection.date = as.Date(Collection.date),
                              region = str_split(Location," / ",simplify = T)[,1],
                              country = str_split(Location," / ",simplify = T)[,2])

metadata.filtered <- metadata.filtered %>% mutate(Pango.lineage2 = ifelse(str_detect(Pango.lineage,"BA\\.1"),"BA.1",
                                                                  ifelse(str_detect(Pango.lineage,"AY\\."),"Delta",
                                                                  as.character(Pango.lineage))))
                                                                  
metadata.filtered <- metadata.filtered %>% filter(Collection.date >= date.start, Collection.date <= date.end)
        
Id.analyzed.v <- metadata.filtered %>% filter(str_detect(Pango.lineage,"BA")) %>% pull(Accession.ID)
                     
#mutation data
mut.interest <- "Spike_L452"
mut.info <- fread(mut.info.name,header=T,sep="\t",check.names=T)
mut.info <- mut.info %>% filter(Id %in% Id.analyzed.v)
mut.info.interest <- mut.info %>% filter(str_detect(mut,mut.interest))


Id.with_Spike_H69_70del.v <- mut.info %>% filter(mut %in% c("Spike_H69del","Spike_V70del")) %>% pull(Id) %>% unique()

metadata.filtered <- metadata.filtered %>% filter(!(Accession.ID %in% Id.with_Spike_H69_70del.v & str_detect(Pango.lineage2,"BA\\.2")))

metadata.filtered <- metadata.filtered %>% left_join(mut.info.interest %>% rename(Accession.ID = Id), by="Accession.ID")
metadata.filtered$mut[is.na(metadata.filtered$mut)] <- ""


pango.interest.v <- c("BA.4","BA.5","BA.2.12.1","BA.2.13","BA.2.9.1","BA.2.11")

metadata.filtered.selected <- metadata.filtered %>% filter(Pango.lineage %in% pango.interest.v)

count.pango.df <- metadata.filtered %>% group_by(Pango.lineage2) %>% summarize(count.pango = n())
count.mut.pango.df <- metadata.filtered %>% filter(mut != "") %>% group_by(Pango.lineage2,mut) %>% summarize(count.mut.pango = n())


count.mut.pango.df <- count.mut.pango.df %>% left_join(count.pango.df,by="Pango.lineage2")
count.mut.pango.df <- count.mut.pango.df %>% mutate(freq.mut.pango = count.mut.pango / count.pango) %>% arrange(desc(freq.mut.pango))

pango.BA2.selected.v <- count.mut.pango.df %>% filter(freq.mut.pango > 0.9, str_detect(Pango.lineage2,"BA\\.2")) %>% pull(Pango.lineage2)

metadata.filtered <- metadata.filtered %>% mutate(Pango.lineage3 = ifelse(!str_detect(Pango.lineage2,"BA\\.2"),as.character(Pango.lineage2),
                                                                   ifelse(Pango.lineage2 %in% pango.BA2.selected.v,as.character(Pango.lineage2),"BA.2")))
                                                                   
country.selected.v <- c("South Africa","USA","Denmark","France","Belgium") #,"Belgium"

metadata.filtered.analyzed <- metadata.filtered %>% filter(country %in% country.selected.v)

count.country.df <- metadata.filtered.analyzed %>% group_by(country) %>% summarize(count.country = n())
count.pango.country.df <- metadata.filtered.analyzed %>% group_by(country,Pango.lineage3) %>% summarize(count.pango.country = n())

count.pango.country.df <- count.pango.country.df %>% inner_join(count.country.df,by="country")
count.pango.country.df <- count.pango.country.df %>% mutate(freq.pango.country = count.pango.country / count.country)

count.pango.country.df <- count.pango.country.df %>% filter(Pango.lineage3 != "Unassigned", count.pango.country >= 100) %>% arrange(country) %>% top_n(5,freq.pango.country) #filter(freq.pango.country > 0.01)
count.pango.country.df %>% as.data.frame()

#min numbers
#limit.prop.recent <- 0.001
limit.count.recent <- 100
#limit.prop.analyzed <- 0.001
limit.count.analyzed <- 100

#Transmissibility
bin.size <- 1
generation_time <- 2.1

#model
multi_nomial_model <- cmdstan_model(stan_f.name)

col.v <- brewer.pal(10, "Paired")
names(col.v) <- unique(count.pango.country.df$Pango.lineage3)

##########transmissibility estimation##########

stat.info.sum <- data.frame()
plot_observed.l <- list()
plot.growth.l <- list()
plot.theta.l <- list()


for(region.interest in country.selected.v){

#region.interest <- "USA"

pango.interest.v <- count.pango.country.df %>% filter(country == region.interest) %>% pull(Pango.lineage3)
pango.interest.v <- c(variant.ref,pango.interest.v[pango.interest.v != variant.ref])


metadata.filtered.analyzed.region <- metadata.filtered %>% filter(country == region.interest, Pango.lineage3 %in% pango.interest.v)

######Transmissibility

metadata.filtered.interest <- metadata.filtered.analyzed.region %>% mutate(date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1, date.bin = cut(date.num,seq(0,max(date.num),bin.size)), date.bin.num = as.numeric(date.bin))
metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin))

metadata.filtered.interest.bin <- metadata.filtered.interest %>% group_by(date.bin.num,Pango.lineage3) %>% summarize(count = n()) %>% ungroup()

metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin %>% spread(key=Pango.lineage3,value = count)
metadata.filtered.interest.bin.spread[is.na(metadata.filtered.interest.bin.spread)] <- 0

X <- as.matrix(data.frame(X0 = 1, X1 = metadata.filtered.interest.bin.spread$date.bin.num))

Y <- metadata.filtered.interest.bin.spread %>% select(- date.bin.num) %>% select(all_of(pango.interest.v))

count.group <- apply(Y,2,sum)
count.total <- sum(count.group)
prop.group <- count.group / count.total

Y <- Y %>% as.matrix()

group.df <- data.frame(group_Id = 1:ncol(Y), group = colnames(Y))

Y_sum.v <- apply(Y,1,sum)

data.stan <- list(K = ncol(Y),
                  N = nrow(Y),
                  D = 2,
                  X = X,
                  Y = Y,
                  generation_time = generation_time,
                  bin_size = bin.size,
                  Y_sum = Y_sum.v)


fit.stan <- multi_nomial_model$sample(
    data=data.stan,
    iter_sampling=2000,
    iter_warmup=500,
    seed=1234,
    parallel_chains = 4,
    #adapt_delta = 0.99,
    max_treedepth = 20,
    chains=core.num)



#growth rate
stat.info <- fit.stan$summary("growth_rate") %>% as.data.frame()
stat.info$Pango.lineage3 <- colnames(Y)[2:ncol(Y)]

stat.info.q <- fit.stan$summary("growth_rate", ~quantile(.x, probs = c(0.025,0.975))) %>% as.data.frame() %>% rename(q2.5 = `2.5%`, q97.5 = `97.5%`)
stat.info <- stat.info %>% inner_join(stat.info.q,by="variable")

stat.info <- stat.info %>% mutate(signif = ifelse(q2.5 > 1,'higher','not higher'))
stat.info <- stat.info %>% arrange(mean)
stat.info <- stat.info %>% mutate(Pango.lineage3 = factor(Pango.lineage3,levels=Pango.lineage3), country = region.interest)

stat.info.sum <- rbind(stat.info.sum,stat.info)

#growth rate plot
draw.df.growth_rate <- fit.stan$draws("growth_rate", format = "df") %>% as.data.frame() %>% select(! contains('.'))
colnames(draw.df.growth_rate) <- colnames(Y)[2:ncol(Y)]

draw.df.growth_rate.long <- draw.df.growth_rate %>% gather(key = class, value = value)

draw.df.growth_rate.long <- rbind(data.frame(class="BA.2",value=1),draw.df.growth_rate.long %>% as.data.frame())

g <- ggplot(draw.df.growth_rate.long,aes(x=class,y=value,color=class,fill=class))
g <- g + geom_hline(yintercept=1, linetype="dashed", alpha=0.5)
g <- g + geom_violin(alpha=0.4,scale="width")
g <- g + stat_summary(geom="pointrange",fun = median, fun.min = function(x) quantile(x,0.025), fun.max = function(x) quantile(x,0.975), size=1,fatten =0.2)

g <- g + scale_color_manual(values = col.v, breaks = names(col.v))
g <- g + scale_fill_manual(values = col.v, breaks = names(col.v))


g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8))
g <- g + xlab('') + ylab('Relative effective reproduction number\n(/BA.2)')
g <- g + theme(legend.position = 'none')
g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
g <- g + ggtitle(region.interest)
g <- g + scale_y_continuous(limits=c(0.6,1.4),breaks=c(0.6,0.8,1,1.2,1.4))
g

plot.growth.l[[region.interest]] <- g


#plot observed

g <- ggplot(metadata.filtered.interest,aes(x=Collection.date,fill=Pango.lineage3))
g <- g + geom_bar(stat = 'count', position = "fill") 
g <- g + scale_x_date(date_labels = "%y-%m", date_breaks = "1 months", date_minor_breaks = "1 month")
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(size=8)
    )
g <- g + scale_color_manual(values = col.v, breaks = names(col.v))

g <- g + ggtitle(region.interest)
g

plot_observed.l[[region.interest]] <- g


#theta

data.freq <- metadata.filtered.interest.bin %>% group_by(date.bin.num) %>% mutate(freq = count / sum(count)) %>% rename(date_Id = date.bin.num, group = Pango.lineage3)



draw.df.theta <- fit.stan$draws("theta", format = "df") %>% as.data.frame() %>% select(! contains('.'))
draw.df.theta.long <- draw.df.theta %>% gather(key = class, value = value)
draw.df.theta.long <- draw.df.theta.long %>% mutate(date_Id = str_match(class,'theta\\[([0-9]+),[0-9]+\\]')[,2] %>% as.numeric(),
                                                            group_Id = str_match(class,'theta\\[[0-9]+,([0-9]+)\\]')[,2] %>% as.numeric())
draw.df.theta.long <- draw.df.theta.long %>% inner_join(group.df,by="group_Id")
draw.df.theta.long <- draw.df.theta.long %>% mutate(Collection.date = as.Date(date_Id,origin=date.start))
draw.df.theta.long.sum <- draw.df.theta.long %>% group_by(date_Id,Collection.date,group) %>% summarize(mean = mean(value),ymin = quantile(value,0.025),ymax = quantile(value,0.975))
draw.df.theta.long.sum <- draw.df.theta.long.sum %>% inner_join(data.freq,by=c("date_Id","group"))

g <- ggplot(draw.df.theta.long.sum,aes(x=Collection.date, y = mean, fill=group, color = group))
g <- g + geom_point(aes(y=freq,size=count), alpha =0.4)
g <- g + geom_ribbon(aes(ymin=ymin,ymax=ymax),alpha=0.4, color=NA)
g <- g + geom_line(size=0.3)
g <- g + scale_x_date(date_labels = "%y-%m", date_breaks = "1 months", date_minor_breaks = "1 month")
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(size=8)
    )
g <- g + scale_color_manual(values = col.v, breaks = names(col.v))
g <- g + scale_fill_manual(values = col.v, breaks = names(col.v))
g <- g + scale_size_continuous(range = c(0.2, 3))
g <- g + ggtitle(region.interest)
g

plot.theta.l[[region.interest]] <- g

}





#write output
write.table(stat.info.sum,txt.growth.rate.name,col.names=T,row.names=F,sep="\t",quote=F)

#plot growth rate
pdf(pdf.growth.rate.name,width=12,height=5)
wrap_plots(plot.growth.l,nrow=1)
dev.off()

#plot observed
pdf(pdf.observed.name,width=20,height=3)
wrap_plots(plot_observed.l,nrow=1)
dev.off()

#plot theta
pdf(pdf.theta.name,width=20,height=3)
wrap_plots(plot.theta.l,nrow=1)
dev.off()

