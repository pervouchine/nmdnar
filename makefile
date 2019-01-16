.PHONY : all

all :: data/upf1xrn1/upf1xrn1vscontrol.tsv data/upf1xrn1/smg6xrn1vscontrol.tsv
all :: hub/shRNA-KD.bed hub/eCLIP_peaks.bed hub/upf1xrn1.bed

##### data download ######

data/shRNA_table.tsv data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv :
	wget https://cb.skoltech.ru/~dp/papers/nmd/data.tar.gz
	tar -xf data.tar.gz
	rm -f data.tar.gz

##### annotation ####

# genes.bed is a bed file containing gene names
data/genes.bed : data/gencode.v19.annotation.gtf
	awk '$$3=="gene"' data/gencode.v19.annotation.gtf | grep 'gene_type "protein_coding";' | perl Perl/print_gff_attributes.pl GEN gene_name | sort -u > data/genes.bed 

# exon_gene.tsv is a two-column tsv file containing exon_id and gene_id 
# only protein-coding and NMD transcripts
data/exon_gene.tsv : data/gencode.v19.annotation.gtf
	awk '$$3=="exon"' data/gencode.v19.annotation.gtf | grep 'gene_type "protein_coding";' | grep 'transcript_type "protein_coding"\|transcript_type "nonsense_mediated_decay"' | perl Perl/print_gff_attributes.pl INT gene_name | sort -u > data/exon_gene.tsv 


#### eCLIP ####

# exon_peaks.bed contains eCLIP peaks in 1000-nt window around exons
# cols 1-6 are exon windows (#4 is exon id); #7 is gene name
# cols 7-17 are eCLIP peaks
data/eCLIP/exon_peaks.bed : data/eCLIP/dump.tsv data/exon_gene.tsv
	awk -v OFS="\t" -v w=1000 '{split($$1,a,"_");if(a[2]<w){a[2]=w};print a[1],a[2]-w,a[3]+w,$$1,1,a[4],$$2}' data/exon_gene.tsv | intersectBed -a stdin -b data/eCLIP/dump.tsv -wa -wb -s > data/eCLIP/exon_peaks.bed

# self_peaks.tsv contains three columns: exon_id, peak_id, gene name
data/eCLIP/self_peaks.bed : data/eCLIP/exon_peaks.bed
	intersectBed -a data/eCLIP/dump.tsv -b data/genes.bed -wa -wb -s | awk '{split($$4,a,"_");if(a[1]==$$17){print}}' | sort -k1,1 -k2,2n > data/eCLIP/self_peaks.bed

# this is a bed file for track hub with all self peaks
hub/eCLIP_peaks.bed :  data/eCLIP/self_peaks.bed
	cut -f1-6 data/eCLIP/self_peaks.bed | awk -v OFS="\t" '{split($$4,a,"_");$$4=a[1];print}' | bedtools merge -s -c 4 -o distinct -i stdin | awk -v OFS="\t" '{print $$1,$$2,$$3,$$5,1000,$$4}' | sort -k1,1 -k2,2n | awk 'BEGIN{print "track name=\"eCLIP\" description=\"eCLIP peaks\" visibility=3 itemRgb=\"Off\""}{print}' > hub/eCLIP_peaks.bed

#### shRNA-KD ####

# this is a make file to compute deltaPSI for all shRNA-KD
shRNA.mk : data/shRNA_table.tsv
	awk -v dir=data/shRNA/ '{out=dir"B07/"$$6"_"$$5".tsv "; print  out ": "dir"A07/"$$1".A07.tsv "dir"A07/"$$2".A07.tsv "dir"A07/"$$3".A07.tsv "dir"A07/"$$4".A07.tsv\n\tmkdir -p "dir"B07/\n\tRscript deltaPSIc.r "dir,$$1,$$2,$$3,$$4,$$5,$$6"\n"; all=all out}END{print "all :" all}' data/shRNA_table.tsv > shRNA.mk

crispr.mk : data/crispr_table.tsv
	awk -v dir=data/crispr/ '{out=dir"B07/"$$6"_"$$5".tsv "; print  out ": "dir"A07/"$$1".A07.tsv "dir"A07/"$$2".A07.tsv "dir"A07/"$$3".A07.tsv "dir"A07/"$$4".A07.tsv\n\tmkdir -p "dir"B07/\n\tRscript deltaPSIc.r "dir,$$1,$$2,$$3,$$4,$$5,$$6"\n"; all=all out}END{print "all :" all}' data/crispr_table.tsv > crispr.mk

# this file pools all deltaPSI for exons in RBPs
data/shRNA/deltaPSI.tsv : shRNA.mk
	make -f shRNA.mk all
	awk '$$1==$$4' data/shRNA/B07/*.tsv > data/shRNA/deltaPSI.tsv

# this is a bed file for track hub
hub/shRNA-KD.bed : data/shRNA/deltaPSI.tsv
	mkdir -p hub/
	awk '$$6>0.05 || $$6<-0.05' data/shRNA/deltaPSI.tsv | awk -v OFS="\t" '{split($$3,a,"_");print a[1],a[2],a[3],"dPSI("$$1","$$2")="$$6,1000,a[4],a[2],a[3], ($$6>0? "255,0,0" : "0,0,255")}' | sort -k1,1 -k2,2n | awk 'BEGIN{print "track name=\"shRNA-KD\" description=\"delta PSI in RBP shRNA-KD\" visibility=3 itemRgb=\"On\""}{print}'> hub/shRNA-KD.bed

####### UPF1/XRN1 KD ######

# this step generates a pdf, png, and tsv table for UPF1/XRN1 KD vs control
data/upf1xrn1/upf1xrn1vscontrol.tsv : data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275416Aligned.sortedByCoord.out.A07.tsv plot.r data/exon_gene.tsv
	Rscript plot.r data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275416Aligned.sortedByCoord.out.A07.tsv 0.001 UPF1/XRN1 data/upf1xrn1/upf1xrn1vscontrol 

# same for SMG6
data/upf1xrn1/smg6xrn1vscontrol.tsv : data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275415Aligned.sortedByCoord.out.A07.tsv plot.r data/exon_gene.tsv
	Rscript plot.r data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275415Aligned.sortedByCoord.out.A07.tsv 0.001 SMG6/XRN1 data/upf1xrn1/smg6xrn1vscontrol

# this is a bed file for track hub
hub/upf1xrn1.bed : data/upf1xrn1/upf1xrn1vscontrol.tsv
	tail -n+2 data/upf1xrn1/upf1xrn1vscontrol.tsv  | awk '$$4>0.05 || $$4<-0.05' | awk -v OFS="\t" '{split($$1,a,"_");print a[1],a[2],a[3],"dPSI(UPF1&XRN1)="$$4,1000,a[4],a[2],a[3], ($$4>0? "255,0,0" : "0,0,255")}' | sort -k1,1 -k2,2n | awk 'BEGIN{print "track name=\"UPF1/XRN1-KD\" description=\"delta PSI in UPF1-KD\" visibility=3 itemRgb=\"On\""}{print}'> hub/upf1xrn1.bed

