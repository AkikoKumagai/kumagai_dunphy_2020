Figure 4
================

### Figure 4A  

``` r
library(dplyr)
library(ggpubr)
library(ggsci)
library(cowplot)

MTBP <- read.delim("../MTBP_peaks.txt", header = TRUE)

# Add distance to Enhancer data.
DistancetoEnhancer <- read.delim("MTBPtoEnhancer.bed", header = F) %>% 
    select(V4, V8) %>% rename(PeakID = V4, Distance.to.Enhancer = V8) %>%
    distinct(PeakID, .keep_all = TRUE)
MTBP4 <- MTBP %>% left_join(DistancetoEnhancer, by = "PeakID") 

# Add information on G4.
G4H <- read.delim("MTBP_G4H.bed", header = FALSE) %>% select(PeakID = V4)
G4H$G4 <- "G4"
MTBP5 <- MTBP4 %>% left_join(G4H, by = "PeakID")

MTBP5$G4[is.na(MTBP5$G4)] <- "no"

# Add information on AP-1 to the MTBP file.
AP1 <- read.delim("MTBP_AP1.bed", header = F) %>% select(PeakID = V4)
AP1$AP1 <- "TGAGTCA"

MTBP6 <- MTBP5 %>% left_join(AP1, by = "PeakID")
MTBP6$AP1[is.na(MTBP6$AP1)] <- "no"

df_TSS_10000 <- MTBP6 %>% filter(Distance.to.TSS > -10000 & Distance.to.TSS < 10000)

df_Enhancer_20000 <- MTBP6 %>% 
    filter(Distance.to.Enhancer > -20000 & Distance.to.Enhancer < 20000)

a <-ggdensity(df_TSS_10000, x = "Distance.to.TSS", y = "..density..", fill = "G4",
              color ="G4", alpha = 0.5, xlim = c(-10000, 10000),
              palette = "startrek", xlab = "Distance to TSS (bp)") 
a2 <- ggpar(a, font.x = c(12, "bold"), font.y = c(14, "bold"), 
            font.xtickslab = c(12, "bold"),  font.ytickslab = c(12, "bold"),
            ylab = "Density", font.legend = c(14, "bold"))
b <- ggdensity(df_Enhancer_20000, x = "Distance.to.Enhancer", y = "..density..", 
               fill = "G4", color ="G4", alpha = 0.5, palette = "startrek", 
               xlim = c(-20000, 20000), xlab = "Distance to Enhancer (bp)")
b2 <- ggpar(b, font.x = c(12, "bold"), font.y = c(14, "bold"), 
            font.xtickslab = c(12, "bold"),  font.ytickslab = c(12, "bold"), 
            ylab = "Density",font.legend = c(14, "bold"))
# Figure 4A top
ggarrange(a2, b2, common.legend = TRUE, legend = "right")

c <- ggdensity(df_TSS_10000, x = "Distance.to.TSS", y = "..density..", 
               fill = "AP1", color ="AP1", alpha = 0.5, palette = "jco", 
               xlim = c(-10000, 10000), xlab = "Distance to TSS (bp)")
d <- ggdensity(df_Enhancer_20000, x = "Distance.to.Enhancer", y = "..density..", 
               fill = "AP1", color ="AP1", alpha = 0.5, palette = "jco", 
               xlim = c(-20000, 20000), xlab = "Distance to Enhancer (bp)")
c2 <- ggpar(c, font.x = c(12, "bold"), font.y = c(14, "bold"), 
            font.xtickslab = c(12, "bold"),  font.ytickslab = c(12, "bold"), 
            ylab = "Density", font.legend = c(14, "bold"))
d2 <- ggpar(d, font.x = c(12, "bold"), font.y = c(14, "bold"), 
            font.xtickslab = c(12, "bold"),  font.ytickslab = c(12, "bold"), 
            ylab = "Density", font.legend = c(14, "bold"))
# Figure 4A bottom
ggarrange(c2, d2, common.legend = TRUE, legend = "right")

# Count the numbers of the samples in each group.
df_TSS_10000 %>% group_by(G4) %>% tally()
df_TSS_10000 %>% group_by(AP1) %>% tally()
df_Enhancer_20000 %>% group_by(G4) %>% tally()
df_Enhancer_20000 %>% group_by(AP1) %>% tally()
```

Perform Kolmogorov-Smirnov test (Supplementary Table S1)

``` r
library(dplyr)
TSS_G4 <- df_TSS_10000 %>% filter(G4 == "G4")
TSS_noG4 <- df_TSS_10000 %>% filter(G4 != "G4")

Enh_G4 <- df_Enhancer_20000 %>% filter(G4 == "G4")
Enh_noG4 <- df_Enhancer_20000 %>% filter(G4 != "G4")

TSS_AP1 <- df_TSS_10000 %>% filter(AP1 == "TGAGTCA")
TSS_noAP1 <- df_TSS_10000 %>% filter(AP1 != "TGAGTCA")

Enh_AP1 <- df_Enhancer_20000 %>% filter(AP1 == "TGAGTCA")
Enh_noAP1 <- df_Enhancer_20000 %>% filter(AP1 != "TGAGTCA")

TSS.G4 <- ks.test(TSS_G4$Distance.to.TSS, TSS_noG4$Distance.to.TSS)
Enh.G4 <- ks.test(Enh_G4$Distance.to.Enhancer, Enh_noG4$Distance.to.Enhancer)

TSS.AP1 <- ks.test(TSS_AP1$Distance.to.TSS, TSS_noAP1$Distance.to.TSS)
Enh.AP1 <- ks.test(Enh_AP1$Distance.to.Enhancer, Enh_noAP1$Distance.to.Enhancer)
```

-----

### Figure 4B  

``` r
# Obtain MTBP summits in Promoter-TSS, Enhancer/Super-enhancer, and Others.
library(readr)
summits <- MTBP6 %>% select(Chr, summit.start, summit.end, Location) 
summits %>% group_by(Location) %>% 
    do(write_tsv(., paste0(unique(.$Location), "_summits.bed"), col_names = FALSE)) 
# This creates three bed files with summits of MTBP peaks.
```

``` bash
computeMatrix reference-point --referencePoint center \
-R Promoter_TSS_summits.bed Enhancer_SuperEnhancer_summits.bed Others_summits.bed \
-S AP1.bw --missingDataAsZero -a 1000 -b 1000 --binSize 10 -out Figure4B.tab.gz

plotProfile -m Figure4B.tab.gz -out Figure4B.pdf --refPointLabel "MTBP" \
--numPlotsPerRow 1 --perGroup --colors darkblue darkblue darkblue
```

-----

### Figure 4C  

``` r
library(tidyr)
t <- MTBP6 %>% group_by(Location, AP1, G4) %>% tally()

all <- MTBP6 %>% group_by(AP1, G4) %>% tally()
all$Location <- "All"

bound <- bind_rows(all, t)

all_data <- all %>% bind_rows(t) %>%  
  mutate(Motifs = case_when(AP1=="no"&G4 =="no" ~ "None", 
                            AP1!="no" & G4!="no" ~ "Both", 
                            AP1!="no"&G4=="no" ~ "AP1", 
                            G4!="no"&AP1=="no"~"G4")) %>% rename(Counts = "n")

all_data$Motifs <- factor(all_data$Motifs,levels = c("G4", "Both", "AP1", "None"))

all_data$Location <- factor(all_data$Location, 
                            levels = c("All", "Promoter_TSS", 
                                       "Enhancer_SuperEnhancer", "Others"))

Fig4C <- ggbarplot(all_data, x="Location", y="Counts",  size = 0, ylab = FALSE,
               xlab=FALSE, order = c("Others", "Enhancer_SuperEnhancer",
                                     "Promoter_TSS", "All"), alpha = 0.8,
               fill = "Motifs", palette = c("#0060FF", "#00b0e0", "#00FFC0",
                                            "#C7c7c7"), color = "white", 
               orientation = "horiz") 
```

-----

### Figure S4C  

``` bash
# Obtain summits in MTBP peaks containing AP-1 motifs. 
# Also obtain matrix of AP-1 motifs near MTBP peaks. 
bedtools intersect -a MTBP_summits.bed -b MTBP_AP1.bed -wa > MTBP_summits_AP1.bed

computeMatrix reference-point --referencePoint center -R MTBP_summits_AP1.bed \
-S AP1.bw --missingDataAsZero -a 1000 -b 1000 --binSize 10 \
-out Matrix_MTBP_AP1.tab.gz 
gunzip Matrix_MTBP_AP1.tab.gz
```

Look at matrix in R.

``` r
library(ggpmisc)
library(tidyr)
library(cowplot)
# Read the matrix file.
TGAGTCA_matrix <- read.delim("../Matrix_MTBP_AP1.tab", skip = 3, header = F)  %>%
    select(-V4, -V5, -V6) %>% 
    rename(Chr = V1, Start = V2, End = V3)

# Separate rows with the AP-1 motif on the right, left, or both sides.
TGAGTCA <- TGAGTCA_matrix %>% 
    mutate(left=select(., V7:V106) %>% apply(1, sum,na.rm=TRUE)) %>%
    mutate(right=select(., V107:V206) %>% apply(1, sum,na.rm=TRUE)) %>% 
    mutate(AP1_site = case_when(left == right ~ "both", left>right ~ "left", 
                                right>left ~ "right")) %>% 
  filter(right != 0 | left !=0)

a <- seq(from = -995, to = 995, by = 10) # Make x-axis 
TGAGTCA_all <- colSums(TGAGTCA[,c(4:203)])/nrow(TGAGTCA)

df <- data.frame(data.frame(a, TGAGTCA_all))  # Combine x-axis and values.

# Plots
FigS4C1 <- ggplot(df, aes(x = a, y = TGAGTCA_all)) + 
    geom_smooth(span = 0.02, se = F, color = "darkblue") + 
    stat_peaks(geom = "text", hjust =1.5, color = "red", ignore_threshold = 0.5) +
    xlab("Distance from MTBP Summits (bp)") + ylab("AP-1 Motif") + 
    theme_cowplot() + panel_border(size = 1, linetype = 1) 

TGAGTCA_group <- TGAGTCA %>% group_by(AP1_site) %>% 
    summarize_at(.vars= c(4:203), list(sum)) 
# Summarize the 4th to 203rd columns to obtain sum. 

TGAGTCA %>% group_by(AP1_site) %>% tally() 
# 862 for both, 3693 for left, and 3663 for right peaks.  10.5 % have both. 89.5% for either side.

TGAGTCA %>% select(c(1:3, 206)) %>% group_by(AP1_site) %>% 
  group_walk(~ write.table(.x, paste0(.y$AP1_site, "_MTBP_TGAGTCA.bed"), 
                           sep="\t", quote=F, col.names = F, row.names=F))  
## right_MTBP_TGAGTCA.bed and left_MTBP_TGAGTCA.bed files were made.

group <- data.frame(t(TGAGTCA_group[, -1]))
colnames(group) <- c("both", "left", "right")

group <- group %>% bind_cols(df) %>% select(1:4) %>% rename(distance = a) 

group_long <- group %>% gather(key = "location", value = "count", -distance) %>%
    mutate(count = count/8218) 
# Transform the data frame to long format. 8218 is the number of rows.

FigS4C2 <- ggplot(group_long, aes(x = distance, y = count, group = location, 
                                  color = location)) + 
    geom_smooth(span = 0.02, se = F) + theme_cowplot() + 
    panel_border( size = 1, linetype = 1) + theme(legend.position=c(0.8, 0.8)) +
    xlab("Distance from MTBP Summits (bp)") + ylab("AP-1 Motif") +
    scale_color_manual(values = c("gray", "red", "blue"))

left <- group_long %>% filter(location == "left")
right <- group_long %>% filter(location  == "right")

FigS4C3 <- ggplot(left, aes(x = distance, y = count, group = location, 
                            color = location)) + 
    geom_smooth(span = 0.02, se = F, color = "red") + 
    theme_cowplot() + panel_border( size = 1, linetype = 1) + ylim(0, 0.020) +
    xlab("Distance from MTBP Summits (bp)") + ylab("AP-1 Motif") 

FigS4C4 <- ggplot(right, aes(x = distance, y = count, group = location, 
                             color = location)) + 
    geom_smooth(span = 0.02, se = F, color = "blue") + 
    theme_cowplot() + panel_border( size = 1, linetype = 1) + 
    ylim(0, 0.020) +
    xlab("Distance from MTBP Summits (bp)") + ylab("AP-1 Motif") 
# bin size = 10  145 bp from the center.
```

-----

### Figure S4D

Look at peaks at TSS to orient peaks.

``` bash
#  Make TSS.bw file.
sort -k1,1 -k2,2n TSS.bed | bedItemOverlapCount hg38 -chromSize=hg38.chrom.sizes \
stdin > TSS.bedgraph

bedGraphToBigWig TSS.bedgraph hg38.chrom.sizes TSS.bw
# Look at TSS placement around Promoter_TSS_summits files.

computeMatrix reference-point --referencePoint center -R Promoter_TSS_summits.bed \
-S TSS.bw --missingDataAsZero -a 1000 -b 1000 --binSize 10 \
-out Matrix_MTBP_TSS.tab.gz 
gunzip Matrix_MTBP_TSS.tab.gz
```

Read the matrix in
R.

``` r
TSS_matrix <- read.delim("../Matrix_MTBP_TSS.tab", skip = 3, header = F)  %>%  
    select(-V4, -V5, -V6) %>% rename(Chr = V1, Start = V2, End = V3)

TSS <- TSS_matrix %>% 
    mutate(left=select(., V7:V106) %>% apply(1, sum,na.rm=TRUE)) %>%
    mutate(right=select(., V107:V206) %>% apply(1, sum,na.rm=TRUE)) %>% 
    mutate(TSS = case_when(left == right ~ "both", left>right ~ "left", 
                           right>left ~ "right")) %>% 
  filter(left != 0 | right !=0)

a <- seq(from = -995, to = 995, by = 10) # Make x-axis 
TSS_all <- colSums(TSS[,c(4:203)])/nrow(TSS)

df <- data.frame(data.frame(a, TSS_all))  # Combine x-axis and values.

FigS4D1 <- ggplot(df, aes(x = a, y = TSS_all)) + 
    geom_smooth(span = 0.05, se = F, color = "darkblue") + 
    stat_peaks(geom = "text", hjust =1.5, color = "red", ignore_threshold = 0.98) +
    xlab("Distance from MTBP Summits (bp)") + ylab("TSS") +
    theme_cowplot() + panel_border(color = "black", size = 1, linetype = 1) 

# Do exactly the same for TSS as for AP-1 motif.
TSS_group <- TSS %>% group_by(TSS) %>% summarize_at(.vars= c(4:203), list(sum)) 
# Summarize the 4th to 203rd columns to get sum. 

TSS %>% group_by(TSS) %>% tally()  # 805 for both, 3331 for left, and 3463 for right peaks.  10.6 % have TSS on both sides. 89.4% on either side.
TSS %>% select(c(1:3, 206)) %>% group_by(TSS) %>% 
    group_walk(~ write.table(.x, paste0(.y$TSS, "_MTBP_TSS.bed"), sep="\t", 
                             quote=F, col.names = F, row.names=F))  

group_TSS <- data.frame(t(TSS_group[, -1]))
colnames(group_TSS) <- c("both", "left", "right")

group_TSS_2 <- group_TSS %>% bind_cols(df) %>% select(1:4) %>% rename(distance = a) 

group_long_TSS <- group_TSS_2 %>% gather(key = "location", value = "count", 
                                         -distance) %>% 
    mutate(count = count/8428) # Transform the data frame to long format.  

FigS4D2 <- ggplot(group_long_TSS, aes(x = distance, y = count, group = location, 
                                      color = location)) + 
    geom_smooth(span = 0.05, se = F) + theme_cowplot() + 
    panel_border(color = "black", size = 1, linetype = 1) +
    theme(legend.position=c(0.8, 0.8)) + xlab("Distance from MTBP Summits (bp)") + 
    ylab("TSS") + scale_color_manual(values = c("gray", "red", "blue"))

left_TSS <- group_long_TSS %>% filter(location == "left")
right_TSS <- group_long_TSS %>% filter(location  == "right")
both_TSS <- group_long_TSS %>% filter(location == "both")

FigS4D3 <- ggplot(left_TSS, aes(x = distance, y = count, group = location,
                                color = location)) + 
    geom_smooth(span = 0.05, se = F, color = "red") + theme_cowplot() +
    panel_border(color = "black", size = 1, linetype = 1) + ylim(0, 0.002) + 
    xlab("Distance from MTBP Summits (bp)") + ylab("TSS") 

FigS4D4 <- ggplot(right_TSS, aes(x = distance, y = count, group = location, 
                                 color = location)) + 
    geom_smooth(span = 0.05, se = F, color = "blue") + theme_cowplot() +
    panel_border(color = "black", size = 1, linetype = 1) + ylim(0, 0.002) +
    xlab("Distance from MTBP Summits (bp)") + ylab("TSS") 
# bin size = 10  145 bp from the center.
```
