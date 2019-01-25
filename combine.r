library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)

df1 = read.delim("data/upf1xrn1/upf1xrn1vscontrol.tsv")
df2 = read.delim("data/upf1xrn1/smg6xrn1vscontrol.tsv")
df3 = merge(df1,df2,by='id',suffixes=c('.UPF1','.SMG6'))
df3$deltaPSIc = with(df3, (deltaPSIc.UPF1+deltaPSIc.SMG6)/2)
df3$p = with(df3, NL.p.UPF1+NL.p.SMG6)

df4 = read.delim("data/shRNA/deltaPSI.tsv", col.names=c('KD','cell','id','gene','deltaPSI','deltaPSIc','z','p','p.adj','KDFC'))[,c('id','deltaPSIc','p','gene')]
df = merge(df3,df4, by='id',suffixes=c('.NMD','.SELF'))
eCLIP= read.delim("data/eCLIP/self_peaks.bed", header=F)
df$self_eCLIP = df$gene %in% eCLIP$V17
df$pval = with(df, p.NMD+p.SELF+10*self_eCLIP)

write.table(df, file="data/combined.tsv", col.names=T, row.names=F, quote=F, sep="\t") 

t = 0.1
b = 0.5

df5 = subset(df, abs(deltaPSIc.SELF)>=t & abs(deltaPSIc.NMD)>=t)
df6 = df5 %>% group_by(gene,deltaPSIc.SELF,deltaPSIc.NMD) %>% slice(which.max(abs(deltaPSIc.SELF*deltaPSIc.NMD)))

p = ggplot(df6, aes(x=deltaPSIc.SELF, y=deltaPSIc.NMD, color=self_eCLIP)) + geom_point(aes(size=sqrt(sqrt(pval)))) 
p = p + geom_label_repel(data=df6, force=10, aes(label = gene), min.segment.length = 1)
p = p + ylim(c(-b,b)) + xlim(c(-b,b)) + ylab(expression(paste(Delta,Psi,"(NMD KD)"))) + xlab(expression(paste(Delta,Psi,"Self KD")))

p = p + geom_segment(aes(x=0,xend=0,y=-b,yend=b),color='black',lty="dashed") + geom_segment(aes(y=0,yend=0,x=-b,xend=b),color='black',lty="dashed")

p = p + theme_classic() + theme(legend.position="none") + scale_colour_manual(values=c("#000088","#880000","#888888"))

print(p)
