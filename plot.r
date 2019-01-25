require(plyr)
require(dplyr)
require(ggplot2)
require(ggrepel)

args <- commandArgs(trailingOnly = TRUE)
# Usage: Rscript plot.r control_file target_file pvalue.cutoff title output

#args = c('data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv','data/upf1xrn1/A07/SRR1275416Aligned.sortedByCoord.out.A07.tsv','0.01','UPF1/XRN1')

compute_psi <- function(df) {
    df$sj_count = df$inc + 2*df$exc
    df$psi   = df$inc / df$sj_count
    subset(df, sj_count>=20)
}

control = compute_psi(read.delim(args[1], col.names=c('id','inc','exc')))
target  = compute_psi(read.delim(args[2], col.names=c('id','inc','exc')))

gene_names = read.delim(pipe("grep -v orf data/exon_gene.tsv | grep -v \"\\.\""), header=F)
colnames(gene_names) = c('id','gene')
merge(merge(control,target,by='id'), gene_names, by='id') -> df

p.cutoff = as.numeric(args[3])
title = args[4]
output = args[5]

df$deltaPSI = round(with(df, psi.y-psi.x), digits=2)
df$log10sjcount = round(with(df, log10(sj_count.x + sj_count.y)), digits=4)
df$log10FC = round(with(df, log10(sj_count.x/sj_count.y)), digits=4)
model = lm(deltaPSI~log10FC, df)
print(summary(model))
df$deltaPSIc = round(model$residuals, digits=2) + p.cutoff/2

df1 = subset(df, deltaPSI!=0)
df1$bin = cut(df1$log10sjcount, breaks=seq(0,6,0.5))
merge(df1, ddply(df1,.(bin), summarise, mean=mean(deltaPSIc), sd=sd(deltaPSIc)))-> df2

df2$z = with(df2, round((deltaPSIc-mean)/sd, digits=2))
df2$p = pnorm(-abs(df2$z))
df2$p.adj = 1 - (1 - df2$p)^dim(df2)[1]
df3 = subset(df2,p.adj<p.cutoff) %>% group_by(gene) %>% slice(which.min(p.adj))

df2$NL.p = round(-log10(df2$p),2)
df2$NL.p.adj = round(-log10(df2$p.adj),2)

p = ggplot(df2, aes(x=log10sjcount, y=deltaPSIc)) + geom_point(size=1, alpha=0.5,aes(color=(p.adj<p.cutoff))) + geom_density_2d(n=100,h=c(0.2,0.08),colour='#666666')
p = p + geom_abline(slope=0,lty="dashed") + geom_label_repel(size=3, data=df3, aes(label = gene,color=(p.adj<p.cutoff)), min.segment.length = 0, force=10)
p = p + xlab("log(SJ count)") + ylab(expression(paste(Delta,Psi,"'(KD-Control)")))
p = p + theme_classic() + theme(legend.position="none") + scale_color_brewer(palette="Set2") + ggtitle(title)

pdf(width=5,height=6, paste(output, ".pdf", sep=""))
print(p)
ggsave(width=5,height=6,paste(output, ".png", sep=""))
dev.off()

fileds = c('id', 'gene', 'deltaPSI', 'deltaPSIc', 'z', 'NL.p', 'NL.p.adj')
write.table(df2[,fileds], file=paste(output, ".tsv", sep=""),col.names=T, row.names=F, quote=F, sep="\t")
