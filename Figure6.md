Figure 6
================

### Processing of the files from EdU-seq experiments

### Making EdU\_seq.bw files (Figure 6A and Figure S6B)

The three replicates of files were aligned to the human genome in the
same manner as CUT\&RUN samples.

``` bash
# Made bw files using Deeptools bamCoverage as follows:
for file in *.bam; do bamCoverage -b ${file} -o bw_files/${file/%.bam/.bw} -bs 50 \
--normalizeUsing RPGC -p max --ignoreDuplicates --effectiveGenomeSize 2913022398 \
--blackListFileName hg38_blacklist.bed
# This makes three bw files (rep1, rep2, rep3) for each of EdU_seq_WT_rep.bw, EdU_seq_deltaC_rep.bw, and EdU_cont_rep.bw. 

# Comparison between WT and control bigwig files and deltaC and control bigwig files. 
# Control is EdU-seq experiments using G1-arrested cells.  
# The values are log2 ratio (sample/control).  bin = 2000
bigwigCompare -b1 EdU_seq_WT_rep1.bw -b2 EdU_cont_rep1.bw -o WT_cont_rep1.bw \
-bs 2000 --blackListFileName hg38_blacklist.bed
# Do this for each replicate. (Also, for each of EdU_seq_deltaC_rep.bw files.)

# Comparison between replicates. 
# Figure S6B
multiBigwigSummary bins -b WT_cont_rep1.bw WT_cont_rep2.bw WT_cont_rep3.bw \
-o multibigwig_WT.npz;

plotCorrelation -in multibigwig_WT.npz --corMethod pearson --whatToPlot scatterplot \
-o FigS6A.pdf --labels rep1 rep2 rep3;

# Merge the files.  The scores are average of three files.
bigWigMerge WT_cont_rep1.bw WT_cont_rep2_WT.bw WT_cont_rep3.bw EdU-seq.bedgraph;

awk 'OFS = "\t"{print $1, $2, $3, $4/3}' EdU-seq.bedgraph | sort -k 1,1 -k2,2n \ 
> EdU-seq_sorted.bedgraph;  
# Score is log2ratio of WT over control averaged between three replicates.

bedGraphToBigWig EdU-seq_sorted.bedgraph hg38.chrom.sizes  EdU_seq_WT.bw
# Do the above for EdU_seq_deltaC.bw
```

-----

### Use HOMER to call peaks.  

``` bash
makeTagDirectory EdU_seq_WT_rep1/ -sspe EdU_seq_WT_rep1.bam;
makeTagDirectory EdU_cont_rep1/ -sspe EdU_cont_rep1.bam;
# Repeat for other two replicates.

findPeaks EdU_seq_WT_rep1/ -i EdU_cont_rep1 -region -size 6000 -minDist 12000 -F 2 >\
EdU1.txt 
# Repeat this for other two replicates.

# Merge peaks of three replicates.
mergePeaks -d given EdU1.txt EdU2.txt EdU3.txt -prefix EdU

grep -v '^#' EdU_EdU1.txt_EdU2.txt EdU_EdU1.txt_EdU3.txt EdU_EdU2.txt_EdU3.txt \
EdU_EdU1.txt_EdU2.txt_EdU3.txt > Initiation_zones.txt

pos2bed.pl Initiation_zones.txt > Initiation_zones.bed # Initiation zones.
```

-----

### Figure S6C

Check the size of the initiation zones.

``` r
library(dplyr)
library(ggpubr)

# Initiation zones = Origins
Origins <- read.delim("Initiation_zones.bed", header = F, skip = 1) %>% mutate(length = V3-V2) %>%
    rename(OriginID = V4, Chr = V1) %>% mutate(size=length/1000) 
# 2473 very early initiation zones.

mean(Origins$length)  # 57.8 kb (OK-seq : mean size = 30 kb, 6-150 kb)
median(Origins$length) # 48 kb

FigS6C <- gghistogram(Origins$size, add = "median",
                      add.params = list(color = "red"), xlim = c(4, 250), 
                      bins = 50, xlab = "Size (kb)", color = "dodgerblue2", 
                      fill = "dodgerblue1",
                      font.label = list(size = 14, face = "bold"))

## Calculate the length of the origins.
total_length <- sum(Origins$size)  # kb
# 143116282  Effective genome size hg38 2913022398
143116282/2913022398 * 100  # 4.9% of the genome.
```

-----

### Figure 6C  

``` bash
# Compare the EdU-seq-WT, EdU-seq-deltaC, and EdU-control read counts. 
multiBigwigSummary BED-file -b EdU_seq_WT_rep1.bw EdU_seq_WT_rep2.bw \ EdU_seq_WT_rep3.bw EdU_seq_deltaC_rep1.bw EdU_seq_deltaC_rep2.bw \ EdU_seq_deltaC_rep3.bw EdU_cont_rep1.bw EdU_cont_rep2.bw EdU_cont_rep3.bw \ --outRawCounts WT_deltaC_EdU.tab -p max --BED Initiation_zones.bed -o EdU_WT_deltaC.npz
# This yields a file (WT_deltaC_EdU.tab) containing log2 average reads for each of the initiation zones.
```

-----

Put the data in R.

``` r
library(tidyverse)
library(ggpubr)
library(ggpmisc)

EdU_WT_deltaC <- read.delim("WT_deltaC_EdU.tab", header = TRUE ) 

EdU_data <- setNames(EdU_WT_deltaC, c("chr", "start", "end", "WT1", "WT2", "WT3",
                                      "deltaC1","deltaC2", "deltaC3", "cont1", 
                                      "cont2", "cont3")) %>% 
    mutate(size = end-start) %>% 
    mutate(WT = (WT1+WT2+WT3)/3, deltaC = (deltaC1+deltaC2+deltaC3)/3, 
           cont = (cont1+cont2+cont3)/3) %>% 
    mutate(WT_cont = WT - cont, deltaC_cont = deltaC-cont) %>% 
    select(chr, start, end, size, WT_cont, deltaC_cont) %>% 
    mutate(ratio = deltaC_cont/WT_cont) %>% mutate(ID = row_number())

mean(EdU_data$WT_cont)
mean(EdU_data$deltaC_cont)

EdU <-gather(EdU_data, group, average_reads, WT_cont:deltaC_cont, factor_key = TRUE)

Fig6C <- ggboxplot(EdU, x = "group", y = "average_reads", ylab = "Average Reads", 
                    fill = "#CD534CFF", outlier.shape = NA, ylim = c(0,7)) +
    scale_y_continuous(expand = c(0, 0)) + 
    stat_compare_means(paired = TRUE, label.y = 5, label.x = 2, label ="p.format")
```

-----

### Figure 6D

Make control bed file (10,000 data points) with regioneR.

``` r
library(regioneR)
control_summits<- createRandomRegions(nregions=10000, length.mean = 2, 
                                      length.sd = 0, genome= "hg38", mask=NULL)
write.table(toDataframe(control_summits), file = "control_summits.bed", 
            sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
```

``` bash
# Initiation zones around MTBP
computeMatrix reference-point --referencePoint center -R Promoter_TSS_summits.bed \
Enhancer_SuperEnhancer_summits.bed Others_summits.bed control_summits.bed  \
-S EdU_seq.bw --missingDataAsZero -a 100000 -b 100000 --binSize 1000 \
-out RO_MTBP2.tab.gz;

plotProfile -m RO_MTBP2.tab.gz -out Fig6D.pdf --regionsLabel "Promoter_TSS" \
"Enhancer_SuperEnhancer" "Others" "random" --samplesLabel EdU-seq --yMin 0 \
--yMax 0.6 --colors green red darkblue yellow --legendLocation upper-left \
--refPointLabel "Center"
```

-----

### Figures 6E and 6F

Look at MTBP peaks around TSS and group them according to the location
of summits of MTBP relative to TSS. (5′side of TSS or 3′ side of TSS)
TSS\_2000.bed file is 1 kb upstream and 1 kb downstream of TSS. MTBP
summits that are within this region are
selected.

``` bash
bedtools intersect -a TSS_2000.bed -b MTBP_summits.bed -wa -wb > TSS_MTBP_peaks.bed
bedtools intersect -a TSS_2000.bed -b MTBP_summits.bed -wa -v > TSS_noMTBP.bed
```

Separate the MTBP peaks with summits on the 5′ side of TSS and 3′ side
of TSS.

``` r
library(dplyr)
library(readr)
TSS_MTBP <- read.delim("../TSS_MTBP_peaks.bed", header = F)

TSS_MTBP_plus <- TSS_MTBP %>% filter(V6 == "+") %>% 
    mutate(loc = case_when(V8 > V2 + 1000 ~ "right", TRUE ~ "left"))
# V8 is the summit of MTBP. "V2 + 1000" is TSS. 

TSS_MTBP_minus <- TSS_MTBP %>% filter(V6 == "-") %>% 
    mutate(loc = case_when(V8 > V3 - 1000 ~ "left", TRUE ~ "right"))

TSS_MTBP_loc <- bind_rows(TSS_MTBP_plus, TSS_MTBP_minus) %>%
    mutate(start = as.integer(V2 + 1000), end = as.integer(start + 1)) %>% 
    select(V1, start, end, V4, V5, V6, loc)

TSS_MTBP_loc_right <- TSS_MTBP_loc %>% filter(loc == "right") 
TSS_MTBP_loc_left <- TSS_MTBP_loc %>% filter(loc == "left")

write.table(TSS_MTBP_loc_right, "TSS_MTBP_loc_right.bed", sep = "\t", 
            col.names = F, row.names = F, quote = F) 
# TSS with MTBP summits on the 3' end.

write.table(TSS_MTBP_loc_left, "TSS_MTBP_loc_left.bed", sep = "\t", 
            col.names = F, row.names = F, quote = F) 
# TSS with MTBP summits on the 5' end.

tss <- read.delim("../TSS.bed", header = F)
```

Draw the graph of MTBP peaks (6E) and EdU-seq reads (6F).

``` bash
# MTBP signal around TSS (Figure 6E)
computeMatrix reference-point --referencePoint TSS -R TSS_MTBP_loc_right.bed \
TSS_MTBP_loc_left.bed -S /WT.bw --missingDataAsZero -a 1000 -b 1000 --binSize 10 \
-out Matrix_TSS_MTBP_rl.tab.gz;

plotProfile -m Matrix_TSS_MTBP_rl.tab.gz -out Figure6D.pdf --regionsLabel "3'" "5'" \
--samplesLabel MTBP --colors darkblue cyan

# Initiation zones around TSS (Figure 6F)
computeMatrix reference-point --referencePoint TSS -p max -R TSS_MTBP_loc_right.bed \
TSS_MTBP_loc_left.bed TSS_noMTBP.bed -S EdU_seq.bw --missingDataAsZero -a 50000 \
-b 50000 --binSize 200 -out TSS_RO_MTBP_rl.tab.gz;

plotProfile -m TSS_RO_MTBP_rl.tab.gz -out Figure6F.pdf \
--regionsLabel "3′" "5′" "no MTBP" --samplesLabel EdU-seq \
--outFileNameData Figure6E.tab
```

Calculate the EdU peak summits from TSS. Sharp peaks around TSS (\<500
bp) and broad peaks around 10 kb are observed.

``` r
library(dplyr)
library(ggpubr)
library(ggpmisc)
library(tidyr)
library(cowplot)
edu_peaks <- read.delim("../Figure6E.tab", header = F, skip = 1)
edu_df <- t(edu_peaks) 
edu <- data.frame(edu_df[c(-1, -2),], stringsAsFactors = FALSE)
edu$X1 <- as.integer(edu_df[3:502,1]) 
edu <- edu %>% mutate(bins = X1-250)
colnames(edu) <- c("X1", "right", "left", "noMTBP", "bins")
edu$right <- as.numeric(edu$right)
edu$left <-as.numeric(edu$left)
edu$noMTBP <- as.numeric(edu$noMTBP)

ggplot(edu, aes(x = bins, y = right)) + 
  geom_smooth(span = 0.01, se = F, color = "darkblue") + 
  stat_peaks(geom = "text", hjust =1.5, color = "red", ignore_threshold = 0.708) + 
  stat_valleys(geom = "text", hjust =1, color = "green", ignore_threshold = 0.3) + 
  xlab("Distance from TSS (bins)") + 
  theme_cowplot() + panel_border(color = "black", size = 1, linetype = 1) 

ggplot(edu, aes(x = bins, y = left)) + 
    geom_smooth(span = 0.01, se = F, color = "darkblue") + 
    stat_peaks(geom = "text", hjust =1, color = "red", ignore_threshold = 0.81) + 
    stat_valleys(geom = "text", hjust =1, color = "green", ignore_threshold = 0.1) + 
    xlab("Distance from TSS (bins)") + theme_cowplot() + 
    panel_border(color = "black", size = 1, linetype = 1) 
# bin = 200
```

-----

### Replication Timing Data

The files for HCT-116 were downloaded from
<https://www2.replicationdomain.com/database.php>. The two files were
merged.

``` bash
bedtools unionbedg -filler ″NA″ -i RT_HCT116_1.bdg RT_HCT116_2.bdg > merged_RT.txt  
# The first few lines were edited out.
awk 'OFS="\t"{print $1, $2, $3, ($4+$5)/2}' merged_RT.txt > RT_HCT116.bdg
```

Segmentation was performed using DNAcopy in Bioconductor.

``` r
library(DNAcopy)

RT <- read.delim("RT_HCT116.bdg", header = F, sep = "\t")
RT$LOC <- as.numeric(RT$V2)
repTiming <- CNA(RT$V4, RT$V1, RT$LOC, data.type ="logratio", sampleid = "RT")

Seg.RT <- segment(repTiming, nperm=10000, alpha=0.01, 
                  undo.splits = "sdundo", undo.SD=65, verbose=2)  
par(ask=T,mar=c(3.1,4.1,1,1))
plot(Seg.RT, plot.type="c")
plot(Seg.RT, plot.type="s")
plot(subset(Seg.RT,chromlist="chr2"), pch=19, pt.cols=c("gray","gray"),
     xmaploc=T,ylim=c(-4.5,4.5))

write.table(Seg.RT$output,"Seg_RT", row.names=F, quote=F, sep= "\t")

RT1  <- Seg.RT$output

# Fill the gaps between the segments.
RT <- RT1 %>% mutate(end = loc.end + 2499, start = loc.start - 2500) %>% 
    mutate(size = end - start)
write.table(RT[,c("chrom", "start", "end", "num.mark", "seg.mean", "size")], 
            "RT.bed", row.names=F, col.names=F, quote=F, sep="\t") 
# RT.bed describes the chr, start, end, RT_scores, and the size of the segments.
```

Intersect RT.bed file with MTBP.bed file and count how many MTBP peaks
are intersected with each replication timing segment.

``` bash
bedClip RT.bed hg38.chrom.sizes RT_c.bed  
# Domains are clipped to the size within chromosomes.  

bedtools intersect -a RT_c.bed -b bed_files/MTBP.bed -c > RT_MTBP_count.bed 
```

For each timing segment, if average RT\_score was more than 1.75, the
segment was designated as “Early”. If RT\_score was less than -1.75, it
was designated as “Late”. If the score was in between, it was designated
as “Mid”.

``` r
library(ggsci)

RT_MTBP <- read.table("RT_MTBP_count.bed", header = F, sep = "\t")
colnames(RT_MTBP) <- c("chr", "start", "end", "data.count", "RT_score",
                       "length", "peak_count")

RT_MTBP$peaks_per_100kb <- (RT_MTBP$peak_count/RT_MTBP$length)*100000

RT_MTBP_T <- RT_MTBP %>% 
    mutate(Timing = case_when(RT_score >= 1.75 ~ "Early", 
                              RT_score < 1.75 & RT_score >= -1.75 ~ "Mid",  
                              RT_score < -1.75 ~ "Late"))

length <- RT_MTBP_T %>% group_by(Timing) %>% tally(length, name = "length") 

peak_count <- RT_MTBP_T %>% group_by(Timing) %>% 
    tally(peak_count, name = "peak_count")

RT_MTBP_T$Timing <- factor(RT_MTBP_T$Timing, levels = c("Early", "Mid", "Late"))

# Add color to the segments for viewing in the browser.  
RT_MTBP_C <- RT_MTBP_T %>% 
    mutate(RGB = case_when(Timing =="Early" ~ "30,144,255", 
                           Timing =="Mid" ~ "255,215,0",
                           Timing =="Late" ~ "169,169,169"))

RT_MTBP_C$strand <- "."

write.table(RT_MTBP_C[,c(1:5,11, 2,3,10)], "RT_color.bed", row.names=F, 
            col.names=F, quote=F, sep="\t") 
write.table(RT_MTBP_T[, c(1:3, 9)], "RT_summary.bed", row.names=F, 
            col.names=F, quote=F, sep="\t") 
```

Intersect with replication timing file and initiation
zone.

``` bash
bedtools intersect -a RT_summary.bed -b bed_files/Initiation_zones.bed  -c > \ RT_EdU_count.bed
```

-----

### Figure S6A  

``` bash
# Add replication timing information to MTBP peak file.
bedtools intersect -a bed_files/MTBP.bed -b RT_summary.bed -wa -wb > MTBP_RT.bed 
```

-----

Add data in R and plot.

``` r
MTBP <- read.delim("MTBP_peaks.txt", header = TRUE)
MTBP_RT <- read.delim("MTBP_RT.bed", header = F) %>% select(PeakID = V4, RT = V14)
MTBP_timing <- MTBP %>% left_join(MTBP_RT, by = "PeakID") %>% select(Location, RT) 
MTBP_timing <- MTBP_timing[complete.cases(MTBP_timing),]

MTBP_timing$Location <- factor(MTBP_timing$Location, levels = c("Promoter_TSS",
                                                            "Enhancer_SuperEnhancer",
                                                            "Others") )

MTBP_timing$RT <- factor(MTBP_timing$RT, levels = c("Early", "Mid", "Late") )
a <- MTBP_timing %>% group_by(RT, Location) %>% tally(n())

FigS6A <- ggbarplot(a, fill="Location", color = "white", y="n", x="RT",
                    position = position_dodge(0.8), 
                    palette = c("#008B45","#EE0000","#3B4992")) 
```

### Figure 6B  

``` r
Fig6B <- ggboxplot(RT_MTBP_T, x="Timing", y= "peaks_per_100kb", 
                   xlab = "Replication Timing", ylab = "No. of Peaks/100kb",  
                   fill = "Timing", ylim = c(0,5), palette = get_palette("jco", 3),
                   outlier.shape=NA) + 
    stat_compare_means(method = "anova", label.x = "Mid", label.y =4) 

# Calculate mean values.
RT_MTBP_T %>% group_by(Timing) %>% tally(mean(peaks_per_100kb))
```

### Figure S6D  

``` bash
computeMatrix reference-point --referencePoint TSS -R TSS_MTBP_loc_right.bed \
TSS_MTBP_loc_left.bed -S bw_files/WT.bw bw_files/deltaC.bw Refs/MNase_hg38.bw H3K4me2_DLD1.bw Refs/H3K4me1*.bigWig Refs/H3K4me2*.bigWig Refs/H3K4me3*.bigWig   Refs/H3K27ac.bw Refs/H3K9Ac*.bigWig Refs/YY1*.bigWig \
--missingDataAsZero -a 1000 -b 1000 --binSize 10 \
-out Matrix_TSS_MTBP_rl_H3.tab.gz;

plotHeatmap -m Matrix_TSS_MTBP_rl_H3.tab.gz -out FigureS6D.pdf 
```
