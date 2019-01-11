.PHONY : all

all :: upf1xrn1vscontrol.pdf smg6xrn1vscontrol.pdf

data/shRNA_table.tsv :
	wget https://cb.skoltech.ru/~dp/papers/nmd/data.tar.gz
	tar -xf data.tar.gz
	rm -f data.tar.gz

data/shRNA/SRSF7_K562.tsv: data/shRNA_table.tsv
	mkdir -p data/shRNA/B07/
	awk '{cmd = "Rscript shRNA.r "$$1" "$$2" "$$3" "$$4" "$$6" data/shRNA/B07/" $$6"_"$$5; system(cmd)}' data/shRNA_table.tsv

#data/crispr/SRSF7_K562.tsv: data/shRNA_table.tsv
#	mkdir -p data/crispr/B07/
#	awk '{cmd = "Rscript shRNA.r "$$1" "$$2" "$$3" "$$4" "$$6" data/shRNA/B07/" $$6"_"$$5; system(cmd)}' data/crispr_table.tsv

all :: data/shRNA/SRSF7_K562.tsv data/crispr/SRSF7_K562.tsv

##############################################

upf1xrn1vscontrol.pdf : data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275416Aligned.sortedByCoord.out.A07.tsv plot.r
	Rscript plot.r data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275416Aligned.sortedByCoord.out.A07.tsv 0.001 UPF1/XRN1 upf1xrn1vscontrol 

smg6xrn1vscontrol.pdf : data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275415Aligned.sortedByCoord.out.A07.tsv plot.r
	Rscript plot.r data/upf1xrn1/A07/SRR1275413Aligned.sortedByCoord.out.A07.tsv data/upf1xrn1/A07/SRR1275415Aligned.sortedByCoord.out.A07.tsv 0.001 SMG6/XRN1 smg6xrn1vscontrol
