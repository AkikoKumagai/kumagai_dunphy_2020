Figure 7
================

### Figure 7A

Make a control file for initiation zones using regioneR.

``` r
library(regioneR)

initiation_zone <- read.delim("bed_files/Initiation_zones.bed", 
                              header = FALSE, skip=1) %>%
    select(c(1,2,3))

control_zone <- randomizeRegions(initiation_zone, genome= "hg38", 
                                 allow.overlaps = FALSE)

write.table(toDataframe(control_zone), file = "control_initiation_zone.bed", 
            sep="\t", col.names = FALSE, row.names = FALSE, quote=FALSE)
```

Count MTBP peaks and YY1 peaks in initiation zones.

``` bash
# Count the number of MTBP peaks per initiation zone and control zone.
bedtools intersect -a bed_files/Initiation_zones.bed -b bed_files/MTBP.bed -c > Init_zone_MTBP.bed
bedtools intersect -a control_initiation_zone.bed -b bed_files/MTBP.bed  \ 
-c > control_zone_MTBP.bed

# Count the number of YY1 peaks per initiation zone and control zone.
bedtools intersect -a bed_files/Initiation_zones.bed -b YY1.bed -c > Init_zone_YY1.bed
bedtools intersect -a control_initiation_zone.bed -b YY1.bed \ 
-c > control_zone_YY1.bed
```

Draw graphs comparing YY1 and MTBP in initiation zones and control
zones.

``` r
library(dplyr)
library(ggpubr)
library(ggsci)

Ini_MTBP <- read.delim("Init_zone_MTBP.bed", header = F) %>% 
    rename(Chr = V1, Start= V2, End = V3, Count = V7) %>% 
    select(Chr, Start, End, Count)

Ini_YY1 <- read.delim("Init_zone_YY1.bed", header = F) %>% 
    rename(Chr = V1, Start= V2, End = V3, Count = V7) %>% 
    select(Chr, Start, End, Count)

Cont_MTBP <- read.delim("control_zone_MTBP.bed", header = F) %>% 
    rename(Chr = V1, Start= V2, End = V3, Count = V4) 

Cont_YY1 <- read.delim("control_zone_YY1.bed", header = F) %>% 
    rename(Chr = V1, Start= V2, End = V3, Count = V4) 

MTBPcount <- bind_rows(Initiation_zone = Ini_MTBP, 
                       Control = Cont_MTBP, .id = "groups") %>% 
    mutate(length = End-Start) %>% 
    mutate(peaks_per_100kb = (Count/length)*100000) 

YY1count <- bind_rows(Initiation_zone = Ini_YY1, 
                      Control = Cont_YY1, .id = "groups") %>% 
    mutate(length = End-Start) %>% 
    mutate(peaks_per_100kb = (Count/length)*100000)

Fig7Al <- ggboxplot(MTBPcount, x="groups", y= "peaks_per_100kb",  
                    ylab = "No. of Peaks/100kb",  fill = "groups",
                    palette = c("#CD534CFF", "#7aa6DCFF"), outlier.shape=NA, 
                    ylim = c(0, 10))  + 
    stat_compare_means(label.y = 8, label.x = 2,label = "p.format")

Fig7Ar <- ggboxplot(YY1count, x="groups", y= "peaks_per_100kb",  
                    ylab = "No. of Peaks/100kb",  fill = "groups",
                    palette = c("#CD534CFF", "#7aa6DCFF"), outlier.shape=NA, 
                    ylim = c(0, 10))  + 
    stat_compare_means(label.y = 8, label.x = 2,label = "p.format")

# Calculate the means.
MTBPcount %>% group_by(groups) %>% tally(mean(peaks_per_100kb))
YY1count %>% group_by(groups) %>% tally(mean(peaks_per_100kb))
```

-----

### Figures 7B and S7A  

``` r
# Obtain the coordinates for the centers of initiation zones.  
sorted_initiation_zone <- initiation_zone %>% 
    mutate(chr = V1, size = V3-V2, center = V2 + as.integer(size/2)) %>% 
    mutate(end = center + 1) %>% 
    select(chr, center, end, size) %>% 
    arrange(desc(size))
# Large to small initiation zones. 

write.table(sorted_initiation_zone, "Initiation_center_sorted.bed", quote = F, 
            sep = "\t", col.names=F, row.names =F )
# This file was used to create Heatmaps.  
```

For heatmaps in Figure S7A, TSS file that covers 2000 bp downstream and
upstream was used to increase the signal.

-----

### Figures 7C and S7B

YY1 HiChIP data is aligned to hg19 genome. MTBP.bw and EdU\_seq data
were converted to hg19. YY1\_HiCHIP data in HCT116 was downloaded from
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2774000>.

``` r
library(regioneR)
library(dplyr)
library(readr)
YY1_HICHIP <- read.csv("GSM2774000_HCT116_YY1_hichip5kb-results.csv.gz",
                       header = TRUE)

# Select high confidence peaks and filter out interchromosomal interactions.
# Also, filter the PETs connecting adjacent bins.
yy1_loops <- YY1_HICHIP %>% filter(Bayes.mixture.1 > 0.9) %>% 
    filter(chromosome1 == chromosome2) %>% filter(end1 != start2)

yy1_loops2 <- yy1_loops %>% mutate(LoopID = 1:n()) # Make LoopID 
yy1_loops2$LoopID <- sprintf('Loop%i', yy1_loops2$LoopID)

yy1_loops3 <- yy1_loops2 %>% select(chromosome1, start1, end2, PET.Count, LoopID) 
# Bed file with start of loop1 and end of loop2.

write.table(yy1_loops3, "yy1_start1_end2.bed", col.names = F, row.names = F, 
            quote = F, sep = "\t")

# Make control regions of the same number and length using regioneR.
original <- read.delim("EdU_hg19.bed", header = F)
random <- randomizeRegions(original, genome= "hg19", allow.overlaps = FALSE)
write.table(toDataframe(random), file = "random.bed", sep="\t", 
            col.names = FALSE, row.names = FALSE, quote=FALSE)
```

Intersect with initiation
zones.

``` bash
bedtools intersect -a EdU_hg19.bed -b yy1_start1_end2.bed -loj > overlapped_YY1.bed;
bedtools intersect -a random.bed -b yy1_start1_end2.bed \ 
-loj > overlapped_YY1_random.bed
```

Count the sum of PET scores per initiation zone.

``` r
detach(package:dplyr)
library(plyr)
library(dplyr)
library(ggpubr)
RO <- read.delim("../overlapped_YY1.bed", header = F) 
random <- read.delim("../overlapped_YY1_random.bed", header = F)

# Count number of PETs.
RO_count <- RO %>% rename(Count = V7) %>% 
    mutate(Count = case_when(V5 < V2-2500  | V6 > V3+2500  ~ "0", 
                             Count == "." ~ "0", 
                             TRUE ~ as.character(Count)))  
# Remove PETs that stick outside of the initiation zone by more than 2500 bp.
# 2500 bp is half of the bin width (5000 bp). 

random_count <- random %>% rename(Count = V7) %>%  
    mutate(Count = case_when(V5 < V2-2500 | V6 > V3+2500  ~ "0", 
                             Count == "." ~ "0", 
                             TRUE ~ as.character(Count)))

RO_count$Count <- as.numeric(RO_count$Count)
random_count$Count <- as.numeric(random_count$Count)

RO_data <- RO_count %>% group_by(V1, V2, V3) %>% summarize(n = sum(Count))

random_data <- random_count %>% group_by(V1, V2, V3) %>% summarize(n = sum(Count))

# Combine the two data frames.
df <- bind_rows(RO = RO_data, random = random_data, .id = "groups" ) %>% 
    mutate(size = V3-V2) %>% mutate(PET_count_10K= (n/size)*10000) %>% 
    filter(size >= 15000) # Small zones are removed.

df %>% group_by(groups) %>% summarise(mean(PET_count_10K))

Fig7C <- ggviolin(df, x = "groups", y = "PET_count_10K", outlier.shape = NA, 
                   fill = "groups", palette = c("#CD534CFF", "#7aa6DCFF"), add = "none", draw_quantiles = 0.5)  + 
    ylim(c(0, 20)) + 
    stat_compare_means(label.y = 15, label.x = 1.5,label = "p.format")
```
