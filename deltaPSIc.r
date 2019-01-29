args <- commandArgs(trailingOnly = TRUE)
print(args)
# Usage: Rscript plot.r target1 target2 control1 control2 pvalue.cutoff title output

require(plyr)
library(qvalue)

compute_psi <- function(df) {
    df$sj_count = df$inc + 2*df$exc
    df$psi   = df$inc / df$sj_count
    subset(df, sj_count>=20)
}

pool_biorep <- function(df1,df2) {
    merge(df1,df2,by="id", all=F) -> df
    df[is.na(df)] <- 0
    with(df, data.frame(id=id, inc=inc.x+inc.y, exc=exc.x+exc.y))
}

dir = args[1]

target1 = read.delim(paste(dir,"A07/",args[2],".A07.tsv", sep=""), col.names=c('id','inc','exc'))
target2 = read.delim(paste(dir,"A07/",args[3],".A07.tsv", sep=""), col.names=c('id','inc','exc'))
target = compute_psi(pool_biorep(target1,target2))

control1 = read.delim(paste(dir,"A07/",args[4],".A07.tsv", sep=""), col.names=c('id','inc','exc'))
control2 = read.delim(paste(dir,"A07/",args[5],".A07.tsv", sep=""), col.names=c('id','inc','exc'))
control = compute_psi(pool_biorep(control1, control2))

gene_names = read.delim(pipe("grep -v orf data/exon_gene.tsv"), header=F)
colnames(gene_names) = c('id','name')
merge(merge(control,target,by='id'), gene_names, by='id') -> df

df$cell = args[6]
df$KD   = args[7]

ddply(subset(df,name==KD),.(),summarize,inc.x=sum(inc.x),inc.y=sum(inc.y)) -> A 
if(dim(A)[1]==0){A=data.frame(inc.x=0,inc.y=0)}
ddply(df,.(),summarize,inc.x=sum(inc.x),inc.y=sum(inc.y)) -> B
df$KDFC=round((A$inc.y/A$inc.x)/(B$inc.y/B$inc.x),digits=2)

output = paste(dir,"B07/",args[7],"_",args[6],".tsv",sep="")
print(output)

df$deltaPSI = round(with(df, psi.y-psi.x), digits=2)
df$log10sjcount = round(with(df, log10(sj_count.x + sj_count.y)), digits=4)
df$log10FC = round(with(df, log10(sj_count.x/sj_count.y)), digits=4)
model = lm(deltaPSI~log10FC, df)
print(summary(model))
df$deltaPSIc = round(model$residuals, digits=2)

df1 = subset(df, deltaPSI!=0)
df1$bin = cut(df1$log10sjcount, breaks=seq(0,6,0.5))
merge(df1, ddply(df1,.(bin), summarise, mean=mean(deltaPSIc), sd=sd(deltaPSIc)))-> df2
df2$z = with(df2, round((deltaPSIc-mean)/sd, digits=2))
df2$p = pnorm(-abs(df2$z))
df2$q = qvalue_truncp(df2$p)$qvalues
df2$p = round(-log10(df2$p),2)
df2$q = round(-log10(df2$q),2)

fileds = c('KD', 'cell', 'id', 'name', 'deltaPSI', 'deltaPSIc', 'z', 'p', 'q','KDFC')
write.table(df2[,fileds], file=output,col.names=T, row.names=F, quote=F, sep="\t")

