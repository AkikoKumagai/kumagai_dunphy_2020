### Processing of Files Used in Figure 3 and 4 Is Described.
***
### Enhancer/Super-enhancer  
Enhancer/Super-enhancer file is from Hnisz et al. (2013) https://www.sciencedirect.com/science/article/pii/S0092867413012270?via%3Dihub#mmc2
The file is hg19. The bed files were lifted over to hg38. The centers of enhancers were determined. "enhancer_center_sorted.bed"
<br>

```bash
bedtools intersect -a MTBP.bed -b refs/Enhancer_HCT116/hg38_Enhancer.bed -wa -u > MTBP_enhancer.bed

bedtools intersect -a MTBP.bed -b refs/Enhancer_HCT116/hg38_SuperEnhancer.bed -wa -u > MTBP_superenhancer.bed

awk 'OFS="\t"{print $4, "Enhancer"}' MTBP_enhancer.bed > MTBP_enhancer.txt

awk 'OFS="\t"{print $4, "SuperEnhancer"}' MTBP_superenhancer.bed > MTBP_superenhancer.txt

cat MTBP_enhancer.txt MTBP_superenhancer.txt > MTBP_SuperEnhancerEnhancer.txt

# Distance to Enhancers from the MTBP peak summits was calculated.
awk 'OFS="\t"{print $1,int($2+($3-$2)/2), int($2+($3-$2)/2)+1}' hg38_Enhancer.bed |bedtools sort > enhancer_center.bed # Center of the enhancers.  

bedtools closest -a MTBP_summits_sorted.bed -b enhancer_center.bed -D a > MTBPtoEnhancer.bed  # Distance from the MTBP peaks to the center of the enhancers.
```
***
### TSS
TSS information was obtained from gencode.v32.basic.annotation.gtf file from Gencode.
https://www.gencodegenes.org/human/


```bash
# Get file containing genes and remove chrM.  Converted from GTF (1-based) to bed files (0-based). 
awk 'OFS="\t"{if ($3 == "gene"){print $1, $4-1, $5, $14, $6, $7 }}' gencode.v32.basic.annotation.gtf | awk 'OFS="\t" {if ($1!="chrM"){print$0}}' > gencode.genes.bed  # 60572 genes.  

# Genes.gtf to TSS.
awk 'OFS="\t" {if($6 == "-") {$2 = $3 - 1} else {$3 = $2 + 1} print}' gencode.genes.bed |awk -F'"' '$1=$1' OFS="\t" | awk 'OFS="\t" {print $1, $2, $3, $4,$6,$7}' >TSS1.bed

# Make bed file covering 1 kb upstream and 1 kb downstream of TSS.
awk 'OFS="\t" {print $1, $2-1000, $3+1000, $4,$5,$6}' TSS.bed > TSS_2000.bed
# 1kb upstream and 1kb downstream.
```
***
### G4 
G4 structures were predicted using G4Hunter.https://doi.org/10.1093/nar/gkw006 using hg38 and threshold = 2.
The results were exported as a bed file. G4H_hg38_2_ref.bed


```bash
# TSS with or without G4. Figure 3C was made using TSS_G4H.bed and TSS_noG4H.bed.
bedtools intersect -a TSS_2000.bed -b G4H_hg38_2_ref.bed -wa -u > TSS_2000_G4H.bed
bedtools intersect -a TSS.bed -b TSS_2000_G4H.bed -wa > TSS_G4H.bed
bedtools intersect -a TSS.bed -b TSS_G4H.bed -wa -v > TSS_noG4H.bed

# MTBP with G4.
bedtools intersect -a bed_files/MTBP.bed -b G4H_hg38_2_ref.bed -wa -u > MTBP_G4H.bed
```
***
### Transcription
RNAseq file (total RNA from DLD1) is from downloaded from GEO (#GSE85688).    
Rokavec M et al. Cancer Res 2017;77(8):1854-1867. PMID: 28130225


```r
library(dplyr)
tss <- read.delim("TSS.bed", header = F) %>% rename(Name =V4)

RNAseq <- read.delim("GSE85685_Rokavec_processed_data.txt", 
                     header = T) %>% 
    select(Name, DLD1) # Column "DLD1" contains the RNA_seq value (RPKM).

tss_RNA <- tss %>% left_join(RNAseq, by = "Name") 
# TSS file and RNAseq file are joined and essential elements were selected.

tss_RNA$rank <- NA
order.rank <- order(tss_RNA$DLD1, decreasing = T)
tss_RNA$rank[order.rank] <- 1:nrow(tss_RNA)  
# Each entry was numbered according to the expression levels.

tss_RNA_select <- tss_RNA[complete.cases(tss_RNA),] %>% arrange(rank)  
#  Rows containing NA were removed.

write.table(tss_RNA_select, "TSS_RNA.bed", quote = F, row.names = F, col.names = F, sep = "\t") # All genes with RNA info were arranged from high to no expression.
```
***
### AP-1 site  
Other motif files (AP1ver2(TTASTCA), RUNX1(AACCACA), TEAD(GGAATGY)) were made in the same manner. AP-1 sites were mapped using HOMER scanMotifGenomeWide.pl


```bash
# Create motif file with zero mismatch. 
seq2profile.pl TGAGTCA 0 AP1 > TGAGTCA.motif  

# Make bed file containing motifs.  
scanMotifGenomeWide.pl TGAGTCA.motif hg38 -bed -mask | awk 'OFS = "\t"{print $1, $2-1, $3,$4,$5,$6}' > AP1.bed 

# Find MTBP peaks with AP1 motif.
bedtools intersect -a bed_files/MTBP.bed -b AP1.bed -wa -u > MTBP_AP1.bed

# Make bigwig file.
sort -k1,1 -k2,2n AP1.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes stdin > AP1.bedgraph;
bedGraphToBigWig AP1.bedgraph hg38.chrom.sizes AP1.bw
```
***
Combine the location information to MTBP_peaks.txt files in R.  

```r
library(dplyr)

MTBP2 <- read.delim("MTBP_peaks.txt", header = TRUE)

Enh <- read.delim("bed_files/MTBP_SuperEnhancerEnhancer.txt", 
                  header = F) %>% rename(PeakID = V1, Enhancer = V2) 

MTBP3 <- MTBP2 %>% left_join(Enh, by="PeakID")  %>% distinct(PeakID, .keep_all = TRUE)

# Change NA to "NO".
MTBP3$Enhancer <- as.character(MTBP3$Enhancer)
MTBP3$Enhancer[is.na(MTBP3$Enhancer)] <- "no"

# Add Location information to the column.
MTBP4 <- MTBP3 %>% 
    mutate(TSS = case_when(Distance.to.TSS < 2000 & Distance.to.TSS > -2000 ~ "TSS",
                           TRUE ~ "no"))  %>%
    mutate(Location = case_when(TSS == "TSS" ~ "Promoter_TSS", 
                                TSS != "TSS" & Enhancer != "no" ~ "Enhancer_SuperEnhancer", TSS != "TSS" & Enhancer == "no" ~ "Others")) %>% 
    select(-TSS, -Enhancer)

write.table(MTBP4, "MTBP_peaks_2.txt", sep = "\t", quote = F, col.names = TRUE,row.names = F )
```
