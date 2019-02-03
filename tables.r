library(xtable)
df1 = read.delim("data/shRNA_table.tsv", header=F)[,-7]
colnames(df1) = c("KD1","KD2","CONTROL1","CONTROL2","cell","target")
df2 = read.delim("data/eCLIP_table.tsv", header=F)
colnames(df2) = c("PEAKS","cell","target","rep")
merge(subset(df2,rep==1), subset(df2,rep==2), by=c("cell","target"),suffixes=c(1,2)) -> df3
merge(df3[,c("cell","target","PEAKS1","PEAKS2")],df1) -> df
tab <- xtable(df, label = 'exptab', caption = 'Accession numbers of eCLIP peaks and KD RNA-seq')
print(tab, tabular.environment = 'longtable', floating = FALSE, include.rownames=F, file="data/exptab.tex")
print(length(unique(df3$target)))

t=0.1
df = read.delim("data/combined.tsv")
arms = c("UPF1","SMG6","NMD","SELF")
df1 = subset(df, p.eCLIP>0 & abs(deltaPSIc.NMD)>=t & abs(deltaPSIc.SELF)>=t & deltaPSIc.NMD*deltaPSIc.SELF<0)
colnames(df1) = c('Exon','Gene',paste("$\\Delta\\Psi_{",arms,"}$",sep=""), paste("p.",c(arms,"eCLIP"), sep=""),"Score")
df1$Exon = gsub("_", ":",as.character(df1$Exon))
tab = xtable(df1, label = 'hits', caption = 'Hits with self-eCLIP and substantial $\\Delta\\Psi$ changes')
print(xtable(tab),include.rownames = FALSE, file="data/hits.tex",type = "latex", sanitize.text.function = function(x){x})


