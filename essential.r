library(ggplot2)
library(ggsignif)

args <- commandArgs(trailingOnly = TRUE)

df = read.delim(args[1])

df = subset(df, df$id %in% read.delim("data/exon_gene.tsv", header=F)$V1)
df = subset(df, df$id %in% read.delim("data/basic_internal_coding_exons.tsv", header=F)$V1)
print(dim(df))

unlist(lapply(strsplit(as.character(df$id),"_"),function(x){as.numeric(x[3])-as.numeric(x[2])+1})) -> df$len
df$phase = df$len %% 3
df$phase=factor(df$phase, levels=c(0,1,2), labels=c('3n','3n+1','3n+2'))
df$essential = factor(df$len %% 3==0, levels=c(T,F), labels=c("Non-essential","Essential"))

p = ggplot(subset(df, deltaPSIc!=0 & abs(deltaPSIc)<0.2), aes(x=phase,fill=essential,y=deltaPSIc)) + geom_boxplot()
p = p + ylim(c(-.2,.25)) + xlab("") + ylab(expression(paste(Delta,Psi,"(KD-Control)"))) + geom_abline(slope=0,lty="dashed")
p = p + geom_signif(comparisons = list(c("3n","3n+1"),c("3n", "3n+2")), map_signif_level=TRUE,y_position = c(.22,.25))
p = p + theme_classic() + theme(legend.position="none") + scale_fill_brewer(palette="Set2")

pdf(width=3,height=4.5,args[2])
print(p)
dev.off()

