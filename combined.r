library(ggplot2)
library(ggrepel)
library(plyr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

nmd = args[1]
self = args[2]
eclip = args[3]
out = args[4]

df1 = na.omit(read.delim(nmd))
df2 = read.delim(self, col.names=c('KD','cell','id','gene','deltaPSI','deltaPSIc','z','p','q','KDFC'))[,c('id','gene','deltaPSIc','p','q')]
f<-function(x) {y=max(abs(x))*sign(x[1]);if(length(x)==2 & prod(x)<0){y=NA};y}
df2 = ddply(df2,.(id,gene),summarize, deltaPSIc=f(deltaPSIc), p=mean(p), q=mean(q))
df = merge(df1,df2, by=c('id','gene'),suffixes=c('.NMD','.SELF'))

ds = read.delim(eclip, col.names=c('id','gene','d','p.eCLIP'))
ds$p.eCLIP = round(ds$p.eCLIP, digits=2)
eCLIP = ds %>% group_by(id,gene,p.eCLIP) %>% slice(which.min(d))
df = merge(df, eCLIP, all.x=T)
df$p.eCLIP[is.na(df$p.eCLIP)] = 0
df$q.eCLIP = df$p.eCLIP
df$Score = with(df, q.NMD + q.SELF + q.eCLIP)

arms = c("UPF1","SMG6","NMD","SELF")
fields = c('id','gene',paste("deltaPSIc",arms, sep="."), paste("p",c(arms,"eCLIP"), sep="."), paste("q",c(arms,"eCLIP"), sep="."),"Score")
write.table(df[, fields], file=paste(out,".tsv",sep=""), col.names=T, row.names=F, quote=F, sep="\t")

pdf(width=7,height=7,paste(out,".pdf",sep=""))
t = 0.1
b = 0.5
df3 = subset(df, abs(deltaPSIc.SELF)>=t & abs(deltaPSIc.NMD)>=t)
df4 = df3 %>% group_by(gene,deltaPSIc.SELF,deltaPSIc.NMD) %>% slice(which.max(abs(deltaPSIc.SELF*deltaPSIc.NMD)))
p = ggplot(df4, aes(x=deltaPSIc.SELF, y=deltaPSIc.NMD, color=(q.eCLIP>0))) + geom_point(aes(size=sqrt(Score))) 
p = p + geom_label_repel(aes(label = gene), min.segment.length=0, force=20)
p = p + ylim(c(-b,b)) + xlim(c(-b,b)) + ylab(expression(paste(Delta,Psi,"(NMD KD)"))) + xlab(expression(paste(Delta,Psi,"Self KD")))
p = p + geom_segment(aes(x=0,xend=0,y=-b,yend=b),color='black',lty="dashed") + geom_segment(aes(y=0,yend=0,x=-b,xend=b),color='black',lty="dashed")
p = p + theme_classic() + theme(legend.position="none") + scale_colour_manual(values=c("#000088","#880000","#888888"))
print(p)
dev.off()
