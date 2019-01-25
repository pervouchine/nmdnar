library(ggplot2)
library(plyr)
library(dplyr)

exon_gene = read.delim("data/exon_gene.tsv", col.names=c('exon','gene'))
peaks= read.delim("data/eCLIP/self_peaks.bed", header=F)[,c(2,3,17)]
colnames(peaks) = c('beg','end','gene')
merge(exon_gene, peaks) -> df
unlist(lapply(strsplit(as.character(df$exon),"_"),function(x){(as.numeric(x[3])+as.numeric(x[2]))/2})) -> df$mid

df$ld = with(df, log10(1+abs(mid - (beg+end)/2)))
df1 = df %>% group_by(gene,exon) %>% slice(which.min(abs(ld)))

merge(df1,ddply(df1,.(gene),summarize,mean=mean(ld),sd=sd(ld))) -> df2
df2$z = with(df2, (ld - mean)/sd)
pdf("data/eCLIP/exon_zscore.pdf")
ggplot(df2,aes(x=gene,y=z)) + geom_boxplot() + coord_flip()
dev.off()

df2$p = round(-log10(pnorm(-abs(df2$z))))
write.table(df2[,c('exon','p')], file="data/eCLIP/exon_pval.tsv", sep="\t", col.names=T, row.names=F, quote=F)
