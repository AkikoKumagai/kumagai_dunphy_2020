Figure 3
================

### Figure 3A  

``` r
# Pie plot for location of MTBP peaks.
library(ggplot2)
library(dplyr)
library(ggsci)
library(scales)
MTBP <- read.delim("MTBP_peaks.txt", header = TRUE)
MTBP %>% group_by(Location) %>% tally()  # The values of this output were used in the next line.

dat <- data.frame(MTBP = c("Promoter-TSS", "Enhancer Super-enhancer", "Others"),
                  value = c(10365, 5257, 14236)) %>% 
    mutate(prop = value/(10365+5257+14236)*100) %>%
    mutate(lab.ypos = cumsum(prop) - 0.5*prop)

dat$prop <- round(dat$prop, digits = 1)

dat$MTBP <- factor(dat$MTBP, levels = c("Others", "Enhancer Super-enhancer",
                                        "Promoter-TSS"))

Fig3A <-ggplot(dat, aes(x="", y = prop, fill = MTBP)) + 
    geom_bar(width = 1, stat = "identity", color = "white") + 
    coord_polar("y", start = 0) + theme_void() +
    geom_text(aes(y = lab.ypos, label = percent(prop/100)), size = 5, 
              color = "white")  + 
    scale_fill_aaas(breaks = c("Promoter-TSS", "Enhancer Super-enhancer", 
                               "Others"))
```

-----

### Figures 3B, 3C, 3D, and S3B

Figures 3B, 3C, and 3D were made using Deeptools computeMatrix and
plotHeatmap. Figure S3B was made using Deeptools computeMatrix and
plotProfile.

``` bash
computeMatrix reference-point --referencePoint center \
-R Promoter_TSS_summits.bed Enhancer_Superenhancer_summits.bed Other_summits.bed \
-S WT.bw deltaC.bw Orc2.bw BG4.bw  --missingDataAsZero -a 3000 -b 3000 --binSize 10 \
-out Matrix_Figure3B.tab.gz

plotHeatmap -m Matrix_Figure3B.tab.gz -out Figure3B.pdf --colorMap YlGnBu YlGnBu \
Blues Blues --refPointLabel "Peak" --sortUsingSamples 4 --whatToShow "heatmap only"
```

Figures 3C and 3D were made using “TSS\_G4H.bed”, “TSS\_noG4H.bed”, and
“TSS\_RNA.bed” files. Figure S3B was made using “TSS.bed”.

### Figure S3C

Correlation between transcription and MTBP peaks at TSS.

``` r
library(dplyr)
library(ggpubr)
MTBP4 <- read.delim("MTBP_peaks.txt", header = T) 
TSS_RNA <- read.delim("TSS_RNA.bed", header = F, sep ="\t") %>% 
    select(c(1:4, 7 ,8)) %>% 
    rename(Gene.Name = V4, Transcription = V7, rank = V8)

MTBP_RNA <- MTBP4 %>% left_join(TSS_RNA, by = "Gene.Name") %>% 
    filter(Location == "Promoter_TSS")  # Filter the MTBP peaks at Promoter_TSS.

data <- MTBP_RNA[complete.cases(MTBP_RNA),]  #  Selected rows with RNA data. 

FigS3C <-ggscatter(data, x = "Transcription", y = "WT.score", color = "Blue", 
                   size = 0.5, alpha = 0.4, ylim = c(80, 1500)) + 
    xscale("log10") + yscale("log10")
```
