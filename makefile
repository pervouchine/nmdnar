.PHONY : all

all :: upf1xrn1vscontrol.pdf smg6xrn1vscontrol.pdf hub/shRNA-KD.bed data/eCLIP/peaks_exon_w1000.bed

data/exon_gene.tsv : data/gencode.v19.annotation.gtf
	awk '$$3=="exon"' data/gencode.v19.annotation.gtf | grep 'gene_type "protein_coding";' | grep 'transcript_type "protein_coding"\|transcript_type "nonsense_mediated_decay"' | print_gff_attributes.pl INT gene_name | sort -u > data/exon_gene.tsv 

data/eCLIP/exon_peaks.bed : data/eCLIP/dump.tsv data/exon_gene.tsv
	awk -v OFS="\t" -v w=1000 '{split($$1,a,"_");if(a[2]<w){a[2]=w};print a[1],a[2]-w,a[3]+w,$$2,1,a[4]}' data/exon_gene.tsv | intersectBed -a stdin -b data/eCLIP/dump.tsv -wa -wb -s > data/eCLIP/exon_peaks.bed

data/eCLIP/peaks_exon_w1000.bed : data/eCLIP/exon_peaks.bed
	awk -v OFS="\t" '{split($$10,a,"_");if(a[1]==$$4){print $$1"_"$$2"_"$$3"_"$$6,a[1]}}' data/eCLIP/exon_peaks.bed | sort -u > data/eCLIP/peaks_exon_w1000.bed

shRNA.mk : data/shRNA_table.tsv
	awk -v dir=data/shRNA/ '{out=dir"B07/"$$6"_"$$5".tsv "; print  out ": "dir"A07/"$$1".A07.tsv "dir"A07/"$$2".A07.tsv "dir"A07/"$$3".A07.tsv "dir"A07/"$$4".A07.tsv\n\tRscript deltaPSIc.r "dir,$$1,$$2,$$3,$$4,$$5,$$6"\n"; all=all out}END{print "all :" all}' data/shRNA_table.tsv > shRNA.mk

data/shRNA/deltaPSI.tsv : shRNA.mk
	make -f shRNA.mk all
	awk '$$1==$$4' data/shRNA/B07/*.tsv > data/shRNA/deltaPSI.tsv

hub/shRNA-KD.bed : data/shRNA/deltaPSI.tsv
	mkdir -p hub/
	awk '$$6>0.05 || $$6<-0.05' data/shRNA/deltaPSI.tsv | awk -v OFS="\t" '{split($$3,a,"_");print a[1],a[2],a[3],"dPSI("$$1")="$$6,1000,a[4],a[2],a[3], ($$6>0? "255,0,0" : "0,0,255")}' | sort -k1,1 -k2,2n | awk 'BEGIN{print "track name=\"shRNA-KD\" description=\"delta PSI in RBP shRNA-KD\" visibility=3 itemRgb=\"On\""}{print}'> hub/shRNA-KD.bed

data/shRNA_table.tsv data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv:
	wget https://cb.skoltech.ru/~dp/papers/nmd/data.tar.gz
	tar -xf data.tar.gz
	rm -f data.tar.gz

##############################################

upf1xrn1vscontrol.pdf : data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275416Aligned.sortedByCoord.out.A07.tsv plot.r data/exon_gene.tsv
	Rscript plot.r data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275416Aligned.sortedByCoord.out.A07.tsv 0.001 UPF1/XRN1 upf1xrn1vscontrol 

smg6xrn1vscontrol.pdf : data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275415Aligned.sortedByCoord.out.A07.tsv plot.r data/exon_gene.tsv
	Rscript plot.r data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275415Aligned.sortedByCoord.out.A07.tsv 0.001 SMG6/XRN1 smg6xrn1vscontrol
