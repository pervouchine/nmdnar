library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

upf1 = args[1]
smg6 = args[2]
self = args[3]
eclip = args[4]
out = args[5]

df1 = read.delim(upf1)
df2 = read.delim(smg6)
df3 = merge(df1,df2,by='id',suffixes=c('.UPF1','.SMG6'))
df3$deltaPSIc = with(df3, (deltaPSIc.UPF1+deltaPSIc.SMG6)/2)
df3$p = with(df3, p.UPF1+p.SMG6)

df4 = read.delim(self, col.names=c('KD','cell','id','gene','deltaPSI','deltaPSIc','z','p','p.adj','KDFC'))[,c('id','deltaPSIc','p','gene')]
df = merge(df3,df4, by='id',suffixes=c('.NMD','.SELF'))

eCLIP= read.delim(eclip, header=F)
df$p.eCLIP = 10*(df$gene %in% eCLIP$V17)
df$pval = with(df, p.NMD + p.SELF + p.eCLIP)

t = 0.1
b = 0.5

df5 = subset(df, abs(deltaPSIc.SELF)>=t & abs(deltaPSIc.NMD)>=t)
df6 = df5 %>% group_by(gene,deltaPSIc.SELF,deltaPSIc.NMD) %>% slice(which.max(abs(deltaPSIc.SELF*deltaPSIc.NMD)))

p = ggplot(df6, aes(x=deltaPSIc.SELF, y=deltaPSIc.NMD, color=(p.eCLIP>0))) + geom_point(aes(size=sqrt(sqrt(pval)))) 
p = p + geom_label_repel(aes(label = gene), min.segment.length=0, force=10)
p = p + ylim(c(-b,b)) + xlim(c(-b,b)) + ylab(expression(paste(Delta,Psi,"(NMD KD)"))) + xlab(expression(paste(Delta,Psi,"Self KD")))

p = p + geom_segment(aes(x=0,xend=0,y=-b,yend=b),color='black',lty="dashed") + geom_segment(aes(y=0,yend=0,x=-b,xend=b),color='black',lty="dashed")

p = p + theme_classic() + theme(legend.position="none") + scale_colour_manual(values=c("#000088","#880000","#888888"))

pdf(width=7,height=7,paste(out,".pdf",sep=""))
print(p)
dev.off()

arms = c("UPF1","SMG6","NMD","SELF")
fields = c('id','gene',paste("deltaPSIc",arms, sep="."),paste("p",c(arms,"eCLIP"), sep="."),"pval")
write.table(df[, fields], file=paste(out,".tsv",sep=""), col.names=T, row.names=F, quote=F, sep="\t")
