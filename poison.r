library(ggplot2)
library(ggsignif)

args <- commandArgs(trailingOnly = TRUE)

data = read.delim(args[1])
exon = read.delim("data/stop.tsv", col.names = c('id','class'))
exon$class = factor(exon$class, levels=c(0,1),labels=c('Control','Poison'))
merge(data, exon, by='id') -> df

p = ggplot(subset(df,abs(deltaPSIc)<.2),aes(x=class,fill=class,y=deltaPSIc)) + geom_boxplot() 
p = p + ylim(c(-.22,.22)) + xlab("") + ylab(expression(paste(Delta,Psi,"(KD-Control)"))) 
p = p + geom_abline(slope=0,lty="dashed") + geom_signif(comparisons = list(c("Poison", "Control")), map_signif_level=TRUE)
p = p + theme_classic() + theme(legend.position="none") + scale_fill_brewer(palette="Set2")

pdf(width=3,height=4.5,args[2])
print(p)
dev.off()
