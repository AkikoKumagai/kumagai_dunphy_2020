---
title: "Figure 5"
output:
  word_document: 
    keep_md: yes
    reference_docx: word-styles-reference-01.docx
  html_document:
    df_print: paged
---

### Figure 5C\

```bash
cat ../AP1.bed ../G4H_hg38_2_ref.bed > AP1_G4.bed
bedtools slop -i ../bed_files/H3K4me2_DLD1.bed -b 200 -g ../../../hdd/refs/hg38.chrom.sizes > H3K4me2_400bp_DLD1.bed

bedtools intersect -a H3K4me2_400bp_DLD1.bed -b AP1_G4.bed -wa -u > H3K4me2_AP1_G4.bed # 33549 peaks for H3K4me2 with G4 or AP-1.
bedtools intersect -a H3K4me2_400bp_DLD1.bed -b AP1_G4.bed -wa -v > H3K4me2_noAP1_noG4.bed # 28955 peaks for H3K4me2 without G4 nor AP1.  
bedtools intersect -a H3K4me2_400bp_DLD1.bed -b ../AP1.bed -wa -u > H3K4me2_AP1.bed;
bedtools intersect -a H3K4me2_400bp_DLD1.bed -b ../G4H_hg38_2_ref.bed -wa -u > H3K4me2_G4.bed;
bedtools intersect -a H3K4me2_400bp_DLD1.bed -b ../bed_files/summits/MTBP_summits.bed -wa -u > H3K4me2_MTBP_summits.bed
```
Put data into R and analyze.

```r
library(tidyverse)
library(ggpubr)

H3K4me2_400 <- read.delim("H3K4me2_400bp_DLD1.bed", header = F) %>% rename(ID = V4)
H3K4me2_AP1 <- read.delim("H3K4me2_AP1.bed", header = F)%>% rename(ID = V4) %>% select(ID)
H3K4me2_AP1$AP1 <- "AP1"

H3K4me2_G4 <- read.delim("H3K4me2_G4.bed", header = F)%>% rename(ID = V4) %>% select(ID)
H3K4me2_G4$G4 <- "G4"

H3K4me2_MTBP_summit <- read.delim("H3K4me2_MTBP_summits.bed", header = F)%>% rename(ID = V4) %>% select(ID)
H3K4me2_MTBP_summit$MTBP <- "MTBP"

H3K4me2 <- H3K4me2_400 %>%  left_join(H3K4me2_AP1, by = "ID") %>%  left_join(H3K4me2_G4, by = "ID")%>%  left_join(H3K4me2_MTBP_summit, by = "ID")
H3K4me2[is.na(H3K4me2)] <- "no"

H3K4me2_motifs <- H3K4me2[, c(1:4, 11:13)] %>% mutate(any_motifs = case_when(AP1=="no" & G4=="no" ~ "without",
                                                                             AP1=="AP1" | G4=="G4" ~ "with"))
# Compare H3K4me2 with or without the motifs for co-localization with MTBP peaks.
H3K4me2_a <- H3K4me2_motifs %>%  group_by(any_motifs, MTBP) %>% tally() %>% rename(Counts = "n")
H3K4me2_a$any_motifs <- factor(H3K4me2_a$any_motifs, levels = c("with", "without"))
H3K4me2_a$MTBP <- factor(H3K4me2_a$MTBP, levels = c("MTBP", "no"))


Fig5C <- ggbarplot(H3K4me2_a, x="any_motifs", y="Counts",  size = 0, ylab = FALSE, xlab = FALSE,
           order = c("with", "without"), alpha = 0.8,
          fill = "MTBP", palette = c("#3333FF", "#009999")) + scale_y_continuous(expand = c(0,0))

mat1<- H3K4me2_a %>% spread(MTBP, Counts) 
mat2 <- mat1[,-1]
row.names(mat2)<- c("with_motif", "without_motif")
# Fisher's exact test
test <- fisher.test(mat2)
test
```
***
### Figure 5D

```bash
computeMatrix reference-point --referencePoint center -p max -R Promoter_TSS_summits.bed  Enhancer_Superenhancer_summits.bed Other_summits.bed -S  ../bw_files/WT.bw ../bw_files/deltaC.bw  ../bw_files/H3K4me2_DLD1.bw  --missingDataAsZero -a 1000 -b 1000 --binSize 10 \
-out Matrix_WT_deltaC_H3K4me2_2.tab.gz

plotHeatmap -m Matrix_WT_deltaC_H3K4me2_2.tab.gz -out Fig5D.pdf --colorMap  YlGnBu YlGnBu Blues Greens --refPointLabel "Peak" --sortUsingSamples 1 --zMax 20 20 80  --whatToShow "heatmap and colorbar"
```
***
### Figures 5E and 5F
Use left_MTBP_TGAGTCA.bed and right_MTBP_TGAGTCA.bed to make a new bed file with strand orientations in column 6.  Make right peaks as "+" and left as "-".

```bash
awk 'OFS = "\t"{print $0,".",".", "+"}' right_MTBP_TGAGTCA.bed > right_MTBP.bed
awk 'OFS = "\t"{print $0,".",".","-"}' left_MTBP_TGAGTCA.bed > left_MTBP.bed

cat right_MTBP.bed left_MTBP.bed | sort -k 1,1 -k2,2n > MTBP_TGAGTCA_dir.bed #7356

awk 'OFS = "\t"{print $0,".",".", "+"}' right_MTBP_TSS.bed > right_TSS_MTBP.bed;
awk 'OFS = "\t"{print $0,".",".","-"}' left_MTBP_TSS.bed > left_TSS_MTBP.bed;

cat right_TSS_MTBP.bed left_TSS_MTBP.bed | sort -k 1,1 -k2,2n > MTBP_TSS_dir.bed
```
***
### Figures S5C and S5D
Use these oriented files (MTBP_TGAGTCA_dir.bed and MTBP_TSS_dir.bed) to make matrix files for Figures S5C and S5D (FigureS5D.tab and FigureS5D.tab)

```bash
# Figure S5C
computeMatrix reference-point --referencePoint center -R MTBP_TGAGTCA_dir.bed \
-S WT.bw MNase_hg38.bw H3K4me2_DLD1.bw H3K4me1.bigWig H3K4me2.bigWig H3K4me3.bigWig H3K9Ac.bigWig \
H3K27ac.bw DNase.bw YY1.bigWig --missingDataAsZero -a 1000 -b 1000 --binSize 10 \
-out direcMTBP.tab.gz

# Figure S5C matrix
plotProfile -m direcMTBP.tab.gz -o dir_MTBP_AP1.pdf --refPointLabel \
"MTBP Peak Summit" --yMin 0 0.06 3 3 3 2 2 0 0     --yMax 12 0.1 40 6 10 8 6 8 1.8 6 \
--regionsLabel " "  --samplesLabel MTBP MNase-seq H3K4me2_DLD1 H3K4me1 H3K4me2 H3K4me3 H3K9ac \
H3K27ac DNase-seq YY1 --outFileNameData FigureS5C.tab 

# Figure S5D
computeMatrix reference-point --referencePoint center -R MTBP_TSS_dir.bed \
-S WT.bw MNase_hg38.bw H3K4me2_DLD1.bw H3K4me1.bigWig H3K4me2.bigWig H3K4me3.bigWig H3K9Ac.bigWig \
H3K27ac.bw DNase.bw YY1.bigWig G4H.bw --missingDataAsZero -a 1000 -b 1000 \
--binSize 10 -out dir_MTBP_TSS.tab.gz

# Figure S5D matrix
plotProfile -m dir_MTBP_TSS.tab.gz -o dir_MTBP_TSS.pdf --yMax 12 0.14 40 4 20 50 20 15 4 15 0.05 --yMin 0 0.05 0 1 9 10 2 5 0 0 0.005  --refPointLabel "MTBP Peak Summit" \ --regionsLabel " " --samplesLabel MTBP MNase-seq H3K4me2_DLD H3K4me1 H3K4me2 H3K4me3 H3K9ac \
H3K27ac DNase-seq YY1 G4 --outFileNameData FigureS5D.tab;
```
***
### Figures 5E and 5F\

```r
library(tidyr)
library(tibble)
library(dplyr)
library(ggpubr)
library(viridis)
df <- read.delim("FigureS5C.tab", skip = 1, header = F) 
data <- df[,c(1,3:202)] 

norm_data <- t(apply(data[2:10,-1], 1, function(x)(x-min(x))/(max(x)-min(x))))  
# The data was normalized from 0 to 1.  

# Transform the matrix to the long format.  
a <- norm_data %>% as_tibble() %>%
  rowid_to_column(var="X") %>%
  gather(key="Y", value="Z", -1) %>%  mutate(Y=as.numeric(gsub("V","",Y)))

Fig5E <- ggplot(a, aes(X, Y, fill= Z)) + 
    geom_tile() +
    theme(legend.position="none") +
    scale_fill_viridis(discrete=FALSE, option = "D", begin = 0, end =0.7) +
    theme_bw()+ coord_flip()+ scale_x_reverse()

# Do the same analysis with TSS.
df_TSS <- read.delim("FigureS5D.tab", skip = 1, header = F) 
data_TSS <- df_TSS[,c(3:202)] 
rownames(data_TSS) <- df_TSS[,1]

# Smooth G4H data.
data_TSS_t <- as.data.frame(t(data_TSS))  # Rows and columns are transposed.

G4H_smooth <- loess(G4 ~ bins, data = data_TSS_t, span = 0.1)
G4H_smoothed <- predict(G4H_smooth)
# Check the smoothing.
ggscatter(data_TSS_t, x = "bins", y = "G4") + geom_smooth(span = 0.1)

# Add smoothed G4 data to the data frame (data_TSS).
newdata <- rbind(data_TSS[1:11,], G4H_smoothed)
rownames(newdata)[12] <- "G4H_smoothed"

# Normalize and draw the graph.
norm_data_TSS <- t(apply(newdata[2:12,], 1, function(x)(x-min(x))/(max(x)-min(x))))  
# The data was normalized from 0 to 1.  
a_TSS <- norm_data_TSS %>% as_tibble() %>%
  rowid_to_column(var="X") %>%
  gather(key="Y", value="Z", -1) %>%  mutate(Y=as.numeric(gsub("V","",Y)))

Fig5F <- ggplot(a_TSS, aes(X, Y, fill= Z)) + 
    geom_tile() +
    theme(legend.position="none") +
    scale_fill_viridis(discrete=FALSE, option = "D", begin = 0, end =0.7) +
    theme_bw()+ coord_flip()+ scale_x_reverse()
```
