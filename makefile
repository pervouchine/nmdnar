.PHONY : all

all :: upf1xrn1vscontrol.pdf smg6xrn1vscontrol.pdf hub/shRNA-KD.bed data/eCLIP/peaks_exon_w1000.bed

data/exon_gene.tsv : data/gencode.v19.annotation.gtf
	awk '$$3=="exon"' data/gencode.v19.annotation.gtf | grep 'gene_type "protein_coding";' | grep 'transcript_type "protein_coding"\|transcript_type "nonsense_mediated_decay"' | print_gff_attributes.pl INT gene_name | sort -u > data/exon_gene.tsv 

data/eCLIP/exon_peaks.bed : data/eCLIP/dump.tsv data/exon_gene.tsv
	awk -v OFS="\t" -v w=1000 '{split($$1,a,"_");if(a[2]<w){a[2]=w};print a[1],a[2]-w,a[3]+w,$$2,1,a[4]}' data/exon_gene.tsv | intersectBed -a stdin -b data/eCLIP/dump.tsv -wa -wb -s > data/eCLIP/exon_peaks.bed

data/eCLIP/peaks_exon_w1000.bed : data/eCLIP/exon_peaks.bed
	awk -v OFS="\t" '{split($$10,a,"_");if(a[1]==$$4){print $$1"_"$$2"_"$$3"_"$$6,a[1]}}' data/eCLIP/exon_peaks.bed | sort -u > data/eCLIP/peaks_exon_w1000.bed

data/shRNA/K562.tsv : data/shRNA/B07/SRSF7_K562.tsv
	awk '$$1==$$4' data/shRNA/B07/*K562.tsv > data/shRNA/K562.tsv

data/shRNA/HepG2.tsv : data/shRNA/B07/SRSF7_HepG2.tsv
	awk '$$1==$$4' data/shRNA/B07/*HepG2.tsv > data/shRNA/HepG2.tsv

hub/shRNA-KD.bed : data/shRNA/K562.tsv data/shRNA/HepG2.tsv
	mkdir -p hub/
	awk '$$6>0.05 || $$6<-0.05' data/shRNA/*.tsv | awk -v OFS="\t" '{split($$3,a,"_");print a[1],a[2],a[3],"dPSI("$$1")="$$6,1000,a[4],a[2],a[3], ($$6>0? "255,0,0" : "0,0,255")}' | sort -k1,1 -k2,2n | awk 'BEGIN{print "track name=\"shRNA-KD\" description=\"delta PSI in RBP shRNA-KD\" visibility=3 itemRgb=\"On\""}{print}'> hub/shRNA-KD.bed

data/shRNA_table.tsv :
	wget https://cb.skoltech.ru/~dp/papers/nmd/data.tar.gz
	tar -xf data.tar.gz
	rm -f data.tar.gz

data/shRNA/B07/SRSF7_K562.tsv: data/shRNA_table.tsv
	mkdir -p data/shRNA/B07/
	awk '{cmd = "Rscript deltaPSIc.r data/shRNA/ "$$1" "$$2" "$$3" "$$4" "$$5" "$$6;system(cmd)}' data/shRNA_table.tsv


all :: data/shRNA/B07/SRSF7_K562.tsv

##############################################

upf1xrn1vscontrol.pdf : data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275416Aligned.sortedByCoord.out.A07.tsv plot.r
	Rscript plot.r data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275416Aligned.sortedByCoord.out.A07.tsv 0.001 UPF1/XRN1 upf1xrn1vscontrol 

smg6xrn1vscontrol.pdf : data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275415Aligned.sortedByCoord.out.A07.tsv plot.r
	Rscript plot.r data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275415Aligned.sortedByCoord.out.A07.tsv 0.001 SMG6/XRN1 smg6xrn1vscontrol
