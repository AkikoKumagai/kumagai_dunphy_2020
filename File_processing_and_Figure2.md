File Processing and Scripts
================
Kumagai and Dunphy (2020)

## File processing

Trim Fastq files.

``` bash
trim_galore --fastqc --paired -q 20 *.fq.gz
```

-----

Align with human index and yeast index.

``` bash
sample=WT_rep1 &&
bowtie2 -p 16 -x GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index \
-1 "${sample}"*val_1.fq.gz -2 "${sample}"*val_2.fq.gz \
-S /sam_files/"${sample}".sam --local --very-sensitive-local --no-unal --no-mixed \
--no-discordant --phred33 -I 10;

bowtie2 -p 16 -x /Bowtie2Index_yeast/genome -1  "${sample}"*val_1.fq.gz \
-2 "${sample}"*val_2.fq.gz -S /spike_in_yeast/"${sample}".sam --local \
--very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 --no-overlap \
--no-dovetail -I 10
```

Remove mtDNA etc. and make bam files and index
file.

``` bash
for file in *.sam; do sed '/chrM/d;/random/d;/chrUn/d' ${file} | samtools view -b \
-q 30  | samtools sort -o /bam_files/${file/%.sam/.bam} 
# Make index.
for file in *.bam; do samtools index ${file} 
```

-----

Count mapped reads.

``` bash
samtools view -F 0x4 WT_rep1.bam | cut -f 1 | sort | uniq | wc -l
```

-----

### BigWig files for MTBP-WT, MTBP-deltaC, and control  

Blacklist file was obtained from
<http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/>.  
Spike-in-ratios (number of reads aligned to human genome)/(number of
reads aligned to yeast genome) were calculated for each experiment.
Scale factors were calculated by using spike-in-ratio (control or
MTBP-ΔC)/spike-in-ratio (MTBP-WT) to normalize the reads to MTBP-WT.
Scale factors ranged between 0.82 and 1.09.  
WT.bw, deltaC.bw, and cont.bw files were made as follows.

``` bash
bamCoverage -b WT_rep1.bam -o WT_rep1.bw -bs 10 --normalizeUsing RPGC \
--ignoreDuplicates --effectiveGenomeSize 2913022398 --scaleFactor 1 \
--blackListFileName  /refs/hg38_blacklist.bed 
## Do the same for each file.  Use own scaleFactors.

# Merge three files for WT, deltaC, and control.
bigWigMerge WT_rep1.bw WT_rep2.bw WT_rep3.bw WT.bedgraph

# Sort and get average scores.
sort -k1,1 -k2,2n WT.bedgraph | \
awk 'BEGIN{OFS = "\t"}($5 = $4/3){print $1,$2,$3,$5}' > WT_adjusted.bedgraph

# Make bw file.
bedGraphToBigWig WT_adjusted.bedgraph /refs/hg38.chrom.sizes WT.bw

# BigWig files for H3K4me2 and control-IgG were made in a similar manner.  
```

-----

### Peak calling with MACS2 for MTBP-WT, MTBP-deltaC, and control  

``` bash
# Call peaks with MACS2 for each replicate.
callpeak -t WT_rep1.bam -c cont_rep1.bam -f BAMPE --nolambda -g 3.05e9 -q 0.01 \
--call-summits -n WT --outdir rep1

callpeak -t deltaC_rep1.bam -c cont_rep1.bam -f BAMPE --nolambda -g 3.05e9 -q 0.01 \
--call-summits -n deltaC --outdir rep1

# Obtain intersection of the replicates. Figure S2.
intervene venn -i rep1/WT_peaks.narrowPeak rep2/WT_peaks.narrowPeak \
rep3/WT_peaks.narrowPeak --save-overlaps -o intervene_three_reps/WT \
--names=rep1, rep2, rep3

intervene venn -i rep1/deltaC_peaks.narrowPeak rep2/deltaC_peaks.narrowPeak \
rep3/deltaC_peaks.narrowPeak --save-overlaps -o intervene_three_reps/deltaC \
--names=rep1, rep2, rep3

# Concatenate the peaks that are in overlapped regions.
# Peaks that are common to at least two out of three replicates are selected.
cat 011_rep2_rep3.bed 101_rep1_rep3.bed 110_rep1_rep2.bed 111_rep1_rep2_rep3.bed \
> WT_rep.bed  

cat 011_rep2_rep3.bed 101_rep1_rep3.bed 110_rep1_rep2.bed 111_rep1_rep2_rep3.bed \
> deltaC_rep.bed

# Call peaks again with three pooled replicates.
callpeak -t WT_rep1.bam WT_rep2.bam WT_rep3.bam -c cont_rep1.bam cont_rep2.bam \
cont_rep3.bam --nolambda -f BAMPE -g 3.05e9 -q 0.01 -n WT --outdir MACS2/WT_pooled

callpeak -t deltaC_rep1.bam deltaC_rep2.bam deltaC_rep3.bam -c cont_rep1.bam \
cont_rep2.bam cont_rep3.bam --nolambda -f BAMPE -g 3.05e9 -q 0.01 -n WT \
--outdir MACS2/deltaC_pooled

# Call peaks for control files.
callpeak -t cont_rep1.bam cont_rep2.bam cont_rep3.bam -c cont_rep1.bam \
cont_rep2.bam cont_rep3.bam --nolambda -f BAMPE -g 3.05e9 -q 0.01 -n WT \
--outdir MACS2/cont_pooled
```

Intersection of pooled peaks and replicated peaks was selected and peaks
in blacklisted area and peaks that intersect with large control peaks
were removed to yield WT.bed and deltaC.bed.

``` bash
bedtools intersect -a WT_pooled/WT_peaks.narrowPeak \
-b intervene_three_reps/WT/sets/WT_rep.bed -wa -u > WT_rep.bed  
# Select all the replicated peaks (29044 peaks) out of 64528 peaks.

bedtools intersect -a WT_rep.bed -b hg38_blacklist.bed -v > WT_noblack.bed 
# 29038 peaks

bedtools intersect -a WT_noblack.bed -b cont_WT.bed -v > WT.bed 
# Do the same for deltaC.
bedtools intersect -a deltaC_noblack.bed -b cont_deltaC.bed -v > deltaC.bed  
```

Peak calling with MACS2 for H3K4me2 and control IgG was performed
without –nolambda. Two replicates for each sample were used and
reproducible peaks between two replicates were selected.

-----

## Figures 2 and S2

### Figure S2B

Selection of reproducible peaks.

``` r
library(ggpubr)
library(dplyr)
library(cowplot)
WT_all <- read.delim("MACS2/WT_pooled/WT_peaks.narrowPeak", header = F)
WT_selected <- read.delim("bed_files/WT.bed", header = F)
WT <- bind_rows("all"=WT_all, "selected"=WT_selected, .id="WT")
Q_6_WT <- WT_selected %>% filter(V9 > 6) # Select peaks that have the -log10qvalue of > 6.
wt_q <- nrow(Q_6_WT)/nrow(WT_selected) # 99.7% of the peaks have -log10qvalue of > 6.
FigS2B_WT <- ggdensity(WT, x = "V9", color = "WT", fill = "WT", palette = c("#00AFBB","#E7B800"), xlim = c(0, 80), ylim = c(0, 5200), y = "..count..", xlab = "-log10qvalue")+
  theme(legend.title = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5)) + ggtitle("WT")

# Do the same for deltaC.
deltaC_all <- read.delim("MACS2/deltaC_pooled/deltaC_peaks.narrowPeak", header = F)
deltaC_selected <- read.delim("bed_files/deltaC.bed", header = F)
deltaC <- bind_rows("all"=deltaC_all, "selected"=deltaC_selected, .id="deltaC")
Q_6_deltaC <- deltaC_selected %>% filter(V9 > 6) # Select peaks that have the -log10qvalue of > 6.
deltaC_q <- nrow(Q_6_deltaC)/nrow(deltaC_selected) # 99.0% of the peaks have -log10qvalue of > 6.
FigS2B_deltaC <- ggdensity(deltaC, x = "V9", color = "deltaC", fill = "deltaC", palette = c("#00AFBB","#E7B800"), xlim = c(0, 80), ylim = c(0, 5200), y = "..count..", 
                           xlab = "-log10qvalue")+
  theme(legend.title = element_blank(), legend.position = c(0.5, 0.8), legend.text = element_text(size = 12), plot.title = element_text(hjust = 0.5)) + ggtitle("deltaC")

plot_grid(FigS2B_WT, FigS2B_deltaC)
wt_q
deltaC_q
```

-----

### Figure 2B

Combine WT.bed and
deltaC.bed.

``` bash
intervene venn -i WT.bed deltaC.bed --save-overlaps -o intervene/WT_deltaC  
# WT_only, WT_deltaC, and deltaC_only files are in intervene/WT_deltaC/sets folder.
```

-----

### Finding peak summits

Peaks were called again with –call-summits in MACS2 as
above.

``` bash
# Intersect final peak files (WT.bed and deltaC.bed) and MACS2 peak files with --call-summits to obtain summits. 
# Many peaks have multiple summits.
bedtools intersect -a bed_files/WT.bed -b MACS2/WT_summits/WT_summits.bed -wa -wb \
> bed_files/WTwithsummits.bed;
bedtools intersect -a bed_files/intervene/WT_deltaC/sets/01_deltaC.bed \
-b MACS2/deltaC_summits/deltaC_summits.bed -wa -wb \
> bed_files/deltaCwithsummits.bed
# deltaC only peaks.
```

Make bed files with peaks containing location information of the
summits.

``` r
# Group by the same ID and take top scored row in the group. 
# Then pick one of the rows for the summits with same value.
WTwithsummits <- read.delim("bed_files/WTwithsummits.bed", header = F) %>% 
    group_by(V4) %>% top_n(1, V15) %>% distinct(V4, .keep_all = TRUE) %>%
    select(V11,V12,V13,V4)

deltaCwithsummits <- read.delim("bed_files/deltaCwithsummits.bed", header = F) %>%
    group_by(V4) %>% top_n(1, V15) %>% distinct(V4, .keep_all = TRUE) %>%
    select(V11,V12,V13,V4)

MTBP_summits <- bind_rows(WTwithsummits, deltaCwithsummits)

write.table(MTBP_summits, "bed_files/MTBP_summits.bed", quote = F, 
            col.names = F, row.names = F, sep = "\t")  
```

Look at the intersection of WT and deltaC. Obtain output of bed files.
These were used to create Figure 2C.

``` bash
# Make bed files that have WT_only, WT_deltaC, and deltaC_only peaks.
bedtools intersect -a WT.bed -b deltaC.bed -u -wa > WT_deltaC.bed;
bedtools intersect -a WT.bed -b deltaC.bed -v > WT_only.bed;
bedtools intersect -a deltaC.bed -b WT.bed -v > deltaC_only.bed; 
bedtools intersect -a WT_summits.bed -b WT_only.bed -wa > WT_only_summits.bed;
bedtools intersect -a WT_summits.bed -b WT_deltaC.bed -wa \
> WT_deltaC_summits.bed;
bedtools intersect -a deltaC_summits.bed -b deltaC_only.bed -wa \
> deltaC_only_summits.bed;
# WT and deltaC
computeMatrix reference-point --referencePoint center -R WT_deltaC_summits.bed \
WT_only_summits.bed deltaC_only_summits.bed \
-S WT.bw deltaC.bw --missingDataAsZero -a 3000 -b 3000 --binSize 10 \
-out Matrix_WT_deltaC.tab.gz

plotHeatmap -m Matrix_WT_deltaC.tab.gz -out hm_WT_deltaC.pdf \
--colorMap YlGnBu --refPointLabel "Peak"  # Figure 2C.
```

Make combined peak file
“MTBP.bed”.

``` bash
cat bed_files/WT.bed bed_files/intervene/WT_deltaC/sets/01_deltaC.bed > \
bed_files/MTBP.bed  ## All peaks.
```

-----

### Figure 2D

Compare MTBP-WT and
MTBP-deltaC.

``` bash
computeMatrix reference-point --referencePoint center -p max -R MTBP_summits.bed \
-S WT.bw deltaC.bw --missingDataAsZero -a 1000 -b 1000 --binSize 10 \
-out Figure2D.tab.gz;

plotProfile -m Figure2D.tab.gz -out Figure2D.pdf  --samplesLabel WT deltaC \
--refPointLabel "MTBP Peak Summits" --colors darkblue darkblue;
```

-----

### Peak Annotation

Annotate peaks was performed using gencode\_basic annotation file
obtained from <https://www.gencodegene>. “MTBP\_peaks\_1.txt” file was
created.

``` bash
annotatePeaks.pl bed_files/MTBP.bed hg38 -gtf gencode.v32.basic.annotation.gtf \
-go ../GO  > MTBP_peaks_1.txt
## Annotate peaks output was also used for Figure S2D.
```

-----

### Figure S2C

Peaks were quantified using Deeptools
multiBigwigSummary.

``` bash
multiBigwigSummary BED-file -b WT.bw cont.bw -out results_WT.npz --BED MTBP.bed \
--outRawCounts WT_counts.tab

multiBigwigSummary BED-file -b deltaC.bw cont.bw -out results_deltaC.npz \
--BED MTBP.bed -p max --outRawCounts deltaC_counts.tab

# Subtract control counts from the samples. Skip header.
tail WT_counts.tab -n+1 | awk 'OFS="\t"{print $1, $2, $3, $4-$5}' > WT-cont.tab
tail deltaC_counts.tab -n+1 | awk 'OFS="\t"{print $1, $2, $3, $4-$5}' > \
deltaC-cont.tab

# Combine the files with MTBP.bed file.
bedtools intersect -a MTBP.bed -b WT-cont.tab -wa -wb > MTBP_count_WT.bed
bedtools intersect -a MTBP.bed -b deltaC-cont.tab -wa -wb > MTBP_count_deltaC.bed
```

In R, calculate the Peak Scores. The output from multiBigwigSummary is
average score in the peak. The bin size = 10.

``` r
library(ggpubr)
library(ggpmisc)

WT_counts <- read.delim("bed_files/MTBP_count_WT.bed", header = F) %>% 
    select(c(1:10, 14)) %>% rename(PeakID = V4, Peak.score = V5) %>% 
    mutate(WT.score = V14*(V3-V2)/10) %>% select(c(1:5), WT.score)

deltaC_counts <- read.delim("bed_files/MTBP_count_deltaC.bed", header = F) %>% 
    select(c(1:10, 14)) %>% rename(PeakID = V4, Peak.score = V5) %>%
    mutate(deltaC.score = V14*(V3-V2)/10) %>% select(PeakID, deltaC.score)

MTBP_counts <- WT_counts %>% left_join(deltaC_counts, by = "PeakID")
MTBP_counts$deltaC.score[MTBP_counts$deltaC.score < 0] <- 0  
# Remove negative values.
MTBP_counts$WT.score[MTBP_counts$WT.score < 0] <- 0

# Regression model with intercept at 0.
WT_deltaC_lm <- lm(deltaC.score ~  0 + WT.score, data = MTBP_counts)
summary(WT_deltaC_lm)

# Sample 5000 data points for the figure.
subset <- sample_n(MTBP_counts, 5000)
FigS2C <- ggscatter(subset, size = 0.5, x = "WT.score", y = "deltaC.score", 
                    color = "dodgerblue", alpha = 0.4, xlim = c(0, 1500), 
                    ylim= c(0, 750)) + 
  geom_abline(slope = WT_deltaC_lm$coefficients, size = 0.2) 
```

-----

### Figure 2E

The difference in read counts per peak between MTBP-WT and MTBP-deltaC.

``` r
library(tidyverse)
library(ggpubr)

MTBP_counts <- read.delim("../MTBP_peaks.txt") %>% select(1, 9, 10) %>% gather(cell_lines, score, 2:3, factor_key = TRUE)

Fig2E <- ggviolin(MTBP_counts, x = "cell_lines", y = "score", ylim = c(0, 1500), draw_quantiles = 0.5, fill = "cell_lines", alpha = 0.5, size =0.5, 
               palette =c("#9BCD9B", "#FFC125"), xlab = "MTBP", width= 0.9, legend = "") +
   stat_summary(fun = "mean", geom="text", aes(label=..y..), position = position_nudge(x=0.5, y=100)) +
    stat_compare_means(label.y =1000, label.x = 0.55, label = "p.format", method = "wilcox", method.args = list(altenative = "greater"), paired = TRUE) 
```

-----

### MTBP\_peak.txt file

“MTBP\_peak.txt” file with summit information and peak scores was
created.

``` r
MTBP1 <- read.delim("MTBP_peaks_1.txt", header = T) %>% select(c(1:4,6,9,10,16))

MTBP1$PeakID <- MTBP1[,1]

# Combine the WT and deltaC peak scores and the locations of the summits.
MTBP1 <- MTBP1 %>% select(PeakID, c(2:9)) 

MTBP_summits <- read.delim("bed_files/summits/MTBP_summits.bed", header = F)  %>%
    rename(PeakID = V4,summit.start = V2,summit.end = V3) %>% select(-V1)

MTBP2 <- MTBP1 %>% 
    left_join(MTBP_counts, by = "PeakID") %>% 
    select(c(1:8, 13, 14)) %>% 
    left_join(MTBP_summits, by = "PeakID") 

write.table(MTBP2, "MTBP_peaks.txt", sep = "\t", quote = F, col.names = TRUE,
            row.names = F )
```
