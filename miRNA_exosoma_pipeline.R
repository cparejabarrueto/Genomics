library(Rhisat2)
#library(heatmap)
library("lattice")
library("tidyr")
library(reshape2)
library(ggplot2)
library(scales)

#wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#gzip -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
#mv Homo_sapiens.GRCh38.dna.primary_assembly.fa genome.fa

tmp <- "/Volumes/MARS4T/"
refs <-  "Homo_sapiens.GRCh38.dna.primary_assembly.fa" #"mature.fa" #"Homo_sapiens.GRCh38.dna.primary_assembly.fa"#"HomoSapiens/ncbi_dataset/rna.fna" #"Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#refs <- "NRG-CP001855.1.fasta"
#refs <- "HS-NC_009800.1.fasta"
hisat2_build(references=refs, outdir=tmp,force=TRUE, prefix="GRCh38")

#hisat2(sequences="1-Exosoma-Negative_S5_L001_R1_001.fastq", index=file.path(tmp, "mature"),type="single", outfile="/Volumes/MARS4T/1-Exosoma-Negative_S5_L001_R1_001.sam", force=TRUE)
#hisat2(sequences="1-Exosoma-Negative_S1_L001_R1_001.fastq", index=file.path(tmp, "mature"),type="single", outfile="/Volumes/MARS4T/1-Exosoma-Negative_S1_L001_R1_001.sam", force=TRUE)
#hisat2(sequences="2-Exosoma-Inf-HS_S7_L001_R1_001.fastq", index=file.path(tmp, "mature"),type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-HS_S7_L001_R1_001.sam", force=TRUE)
#hisat2(sequences="2-Exosoma-Inf-HS_LB_C_S3_L001_R1_001.fastq", index=file.path(tmp, "mature"),type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-HS_LB_C_S3_L001_R1_001.sam", force=TRUE)
#hisat2(sequences="2-Exosoma-Inf-NRG_S6_L001_R1_001.fastq", index=file.path(tmp, "mature"), type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-NRG_S6_L001_R1_001.sam", force=TRUE)
#hisat2(sequences="2-Exosoma-Inf-NRG_LB_B_S2_L001_R1_001.fastq", index=file.path(tmp, "mature"),type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-NRG_LB_B_S2_L001_R1_001.sam", force=TRUE)


#INDEX
hisat2(sequences="1-Exosoma-Negative_S5_L001_R1_00.fastq", index=file.path(tmp, "GRCh38"),type="single", outfile="/Volumes/MARS4T/1-Exosoma-Negative_S5_L001_R1_001.sam", force=TRUE)
hisat2(sequences="1-Exosoma-Negative_S1_L001_R1_00.fastq", index=file.path(tmp, "GRCh38"),type="single", outfile="/Volumes/MARS4T/1-Exosoma-Negative_S1_L001_R1_001.sam", force=TRUE)
hisat2(sequences="2-Exosoma-Inf-HS_S7_L001_R1_00.fastq", index=file.path(tmp, "GRCh38"),type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-HS_S7_L001_R1_001.sam", force=TRUE)
hisat2(sequences="2-Exosoma-Inf-HS_LB_C_S3_L001_R1_00.fastq", index=file.path(tmp, "GRCh38"),type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-HS_LB_C_S3_L001_R1_001.sam", force=TRUE)
hisat2(sequences="2-Exosoma-Inf-NRG_S6_L001_R1_00.fastq", index=file.path(tmp, "GRCh38"), type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-NRG_S6_L001_R1_001.sam", force=TRUE)
hisat2(sequences="2-Exosoma-Inf-NRG_LB_B_S2_L001_R1_00.fastq", index=file.path(tmp, "GRCh38"),type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-NRG_LB_B_S2_L001_R1_001.sam", force=TRUE)

#/Volumes/MARS4T/1-Exosoma-Negative_S5_L001_R1_001.sam /Volumes/MARS4T/1-Exosoma-Negative_S1_L001_R1_001.sam /Volumes/MARS4T/2-Exosoma-Inf-HS_S7_L001_R1_001.sam /Volumes/MARS4T/2-Exosoma-Inf-HS_LB_C_S3_L001_R1_001.sam /Volumes/MARS4T/2-Exosoma-Inf-NRG_S6_L001_R1_001.sam "Volumes/MARS4T/2-Exosoma-Inf-NRG_LB_B_S2_L001_R1_001.sam

#htseq-count --format=sam --order=name --stranded=yes --type=exon --idattr=gene_id --mode=intersection-strict 1-Exosoma-Negative_S5_L001_R1_001.sam HomoSapiens/gencode.v44.tRNAs.gtf > 1.count

# htseq-count -f sam  -s yes  -i gene_id -m intersection-strict /Volumes/MARS4T/1-Exosoma-Negative_S5_L001_R1_001.sam /Volumes/MARS4T/1-Exosoma-Negative_S1_L001_R1_001.sam /Volumes/MARS4T/2-Exosoma-Inf-HS_S7_L001_R1_001.sam /Volumes/MARS4T/2-Exosoma-Inf-HS_LB_C_S3_L001_R1_001.sam /Volumes/MARS4T/2-Exosoma-Inf-NRG_S6_L001_R1_001.sam "Volumes/MARS4T/2-Exosoma-Inf-NRG_LB_B_S2_L001_R1_001.sam HomoSapiens/ncbi_dataset/genomic.gtf > /Volumes/MARS4T/pipeline.count

#hisat2(sequences="1-Exosoma-Negative_S5_L001_R1_001.fastq", index=file.path(tmp, "rna"),type="single", outfile="/Volumes/MARS4T/1-Exosoma-Negative_S5_L001_R1_001.sam", force=TRUE,execute=TRUE)
#hisat2(sequences="1-Exosoma-Negative_S1_L001_R1_001.fastq", index=file.path(tmp, "rna"),type="single", outfile="/Volumes/MARS4T/1-Exosoma-Negative_S1_L001_R1_001.sam", force=TRUE)
#hisat2(sequences="2-Exosoma-Inf-HS_S7_L001_R1_001.fastq", index=file.path(tmp, "rna"),type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-HS_S7_L001_R1_001.sam", force=TRUE)
#hisat2(sequences="2-Exosoma-Inf-HS_LB_C_S3_L001_R1_001.fastq", index=file.path(tmp, "rna"),type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-HS_LB_C_S3_L001_R1_001.sam", force=TRUE)
#hisat2(sequences="2-Exosoma-Inf-NRG_S6_L001_R1_001.fastq", index=file.path(tmp, "rna"),type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-NRG_S6_L001_R1_001.sam", force=TRUE)
#hisat2(sequences="2-Exosoma-Inf-NRG_LB_B_S2_L001_R1_001.fastq", index=file.path(tmp, "rna"),type="single", outfile="/Volumes/MARS4T/2-Exosoma-Inf-NRG_LB_B_S2_L001_R1_001.sam", force=TRUE)

# RUN
#/usr/local/bin/AdapterRemoval --file1 2-Exosoma-Inf-NRG_S6_L001_R1_001.fastq --basename 2-Exosoma-Inf-NRG_S6_L001_R1_00
#.././subread-2.0.3-macOS-x86_64/bin/featureCounts -s 0 -T 2 -t miRNA -g Name -a hsa.gff3 -o counts.txt /Volumes/MARS4T/1-Exosoma-Negative_S5_L001_R1_001.sam /Volumes/MARS4T/1-Exosoma-Negative_S1_L001_R1_001.sam /Volumes/MARS4T/2-Exosoma-Inf-HS_S7_L001_R1_001.sam /Volumes/MARS4T/2-Exosoma-Inf-HS_LB_C_S3_L001_R1_001.sam /Volumes/MARS4T/2-Exosoma-Inf-NRG_S6_L001_R1_001.sam /Volumes/MARS4T/2-Exosoma-Inf-NRG_LB_B_S2_L001_R1_001.sam

test <- read.table(file='counts.txt',header=T, col.names = c("Geneid","Chr","Start","End","Strand","Length","1-Exosoma-Negative_S5_L001_R1_001.sam","1-Exosoma-Negative_S1_L001_R1_001.sam","2-Exosoma-Inf-HS_S7_L001_R1_001.sam","2-Exosoma-Inf-HS_LB_C_S3_L001_R1_001.sam","2-Exosoma-Inf-NRG_S6_L001_R1_001.sam","2-Exosoma-Inf-NRG_LB_B_S2_L001_R1_001.sam"))


colnames(test) <- c("gene","Chr","Start","End","Strand","Length","Negative1","Negative2","HS_LB_SALES","HS_LB","NRG_LB_SALES","NRG_LB")
order1 <- test[order((test$Negative1+test$Negative2+test$HS_LB_SALES+test$HS_LB+test$NRG_LB_SALES+test$NRG_LB),decreasing=TRUE),]
write.table(order1, "order1.txt", sep = "\t", row.names = F, quote = F)
order2 <- read.delim("order1.txt", h = T)

#eliminar 
order2['HS-Negative LB+'] <- order2$HS_LB_SALES-order2$Negative1
order2['NRG-Negative LB+'] <- order2$NRG_LB_SALES-order2$Negative1
order2['Exosome Inf HS LB+']<-scale(order2$`HS-Negative LB+`)
order2['Exosome Inf NRG LB+']<-scale(order2$`NRG-Negative LB+`)

order2['HS-Negative LB'] <- order2$HS_LB-order2$Negative2
order2['NRG-Negative LB'] <- order2$NRG_LB-order2$Negative2
order2['Exosome Inf HS LB']<-scale(order2$`HS-Negative LB`)
order2['Exosome Inf NRG LB']<-scale(order2$`NRG-Negative LB`)

write.table(order2, "order2.txt", sep = "\t", row.names = F, quote = F)
matrix <- data.frame(order2$gene,scale(order2$`HS-Negative LB+`),scale(order2$`NRG-Negative LB+`),scale(order2$`HS-Negative LB`),scale(order2$`NRG-Negative LB`))
colnames(matrix) <- c("gene","Exosome Inf HS LB+","Exosome Inf NRG LB+","Exosome Inf HS LB","Exosome Inf NRG LB")
co=melt(matrix)
head(co)
ggplot(co, aes( variable,gene)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value)) + # background colours are mapped according to the value column
  #geom_text(aes(fill = co$value, label = round(co$value, 2))) + # write the values
  scale_fill_gradient2(low = muted("darkred"), 
                       mid = "white", 
                       high = muted("midnightblue"), 
                       midpoint = 0) + # determine the colour
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 12,face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold")) + 
 # ggtitle("Correlation Plot") + 
  theme(legend.title=element_text(face="bold", size=14)) + 
  scale_x_discrete(name="") +
  scale_y_discrete(name="") +
  labs(fill="Scale")



