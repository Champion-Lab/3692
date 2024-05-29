setwd("~/04_Champion_Lab/02_N-terminal_Acetylation/3740_3692_DDABUP/DDH 3740 All LFQ/")

library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(ggrepel)
library(ggtext)

peptides <- read.csv("db.protein-peptides.csv")

peptides <- peptides[!grepl("CONTAM", peptides$Accession),]
peptides <- peptides[peptides$Start <= 2,]
peptides$Peptide <- str_remove(peptides$Peptide, "^.\\.")
peptides$Peptide <- str_remove(peptides$Peptide, "\\..$")


peptides$NtermAcet_L <- grepl("^[A-Z]\\(\\+42\\.01\\)", peptides$Peptide)
peptides$NtermAcet_H <- grepl("^[A-Z]\\(\\+45\\.03\\)", peptides$Peptide)
peptides$UnLabelled <- !peptides$NtermAcet_H & !peptides$NtermAcet_L

peptidesLong <- peptides  %>% pivot_longer(which(grepl("^Area.x[0-9]", names(peptides))),
                                           values_to = "Area",
                                           names_to = "sample")
drops <- c("USED", "Quality", "Significance","Avg..ppm","X1.K0","Avg..Area",
           "Sample.Profile..Ratio.","Area.WT","Area.Del",             
           "Area.Comp","Group.Profile..Ratio.","Max.Ratio", "X.Vector",
           names(peptides)[grepl("Intensity.", names(peptides))],
           names(peptides)[grepl("X.Spec", names(peptides))])

peptides <- peptides[,!names(peptides) %in% drops]


peptidesLong <- peptides  %>% pivot_longer(which(grepl("^Area.x[0-9]", names(peptides))),
                                           values_to = "Area",
                                           names_to = "sample")


peptidesLong$Inj <- str_remove(peptidesLong$sample, "Area.x[0-9]+L_")
peptidesLong$strain <- str_remove(str_extract(peptidesLong$Inj, "\\_[^_]+$"), "_")
peptidesLong$Biorep <- str_remove(str_extract(peptidesLong$Inj, "[PS][0-9]"), "[PS]")
peptidesLong$sampleType <- str_extract(peptidesLong$Inj, "^[PS]")
peptidesLong$Techrep <- str_remove_all(str_extract(peptidesLong$Inj, "_[0-9]_"), "_")

NacetpeptidesLong <- peptidesLong

NacetpeptidesLong$Accession2 <- str_remove(str_extract(NacetpeptidesLong$Accession, "^[^\\|]*\\|"), "\\|")
NacetpeptidesLong$gene <- str_remove(str_remove(str_extract(NacetpeptidesLong$Accession, "^[^\\|]*\\|[^\\|]*\\|"),
                                                "^[^\\|]*\\|"), "\\|")

accessionsWithNtermAcet <- unique(NacetpeptidesLong$Accession2[NacetpeptidesLong$NtermAcet_L])

NacetpeptidesLong <- NacetpeptidesLong[NacetpeptidesLong$Accession2 %in% accessionsWithNtermAcet,]

NacetpeptidesLong$acetSeq[NacetpeptidesLong$NtermAcet_L] <- str_replace(NacetpeptidesLong$Peptide[NacetpeptidesLong$NtermAcet_L],
                                                                        "\\(\\+42\\.01\\)",
                                                                        "(ac)")

NacetpeptidesLong$acetSeq[NacetpeptidesLong$NtermAcet_H] <- str_replace(NacetpeptidesLong$Peptide[NacetpeptidesLong$NtermAcet_H],
                                                                        "\\(\\+45\\.03\\)",
                                                                        "(ac)")
NacetpeptidesLong$acetSeq[NacetpeptidesLong$UnLabelled] <- NacetpeptidesLong$Peptide[NacetpeptidesLong$UnLabelled]
NacetpeptidesLong$noAcetSeq <- str_remove(NacetpeptidesLong$acetSeq, "\\(ac\\)")

STRAIN <- "3692"
DELSTRAIN <- paste0("d", STRAIN)
COMPSTRAIN <- paste0("c", STRAIN)

sampleType <- "P"

intermediate <- NacetpeptidesLong[NacetpeptidesLong$sampleType == sampleType,]


peps <- intermediate[intermediate$strain == "WT" |
                       intermediate$strain == "cross" |
                            grepl(STRAIN, intermediate$strain),]

peps$category <- ifelse(peps$UnLabelled, "Unlabelled",
                        ifelse(peps$NtermAcet_H, "Heavy", "Light"))

peps2 <- peps %>% group_by(Peptide, category) %>%
  mutate(rows = n(),
         IDs = sum(!is.na(Area)),
         meanArea = mean(Area, na.rm = T))

lights <- peps2[peps2$NtermAcet_L,]

peps3 <- lights %>% group_by(Accession2, category) %>%
  mutate(bestFlierArea = max(meanArea, na.rm = T),
         bestFlier = bestFlierArea == meanArea)

bestFliers <- peps3[peps3$bestFlier &
                      !is.na(peps3$bestFlier),]

bestFliers$HeavyArea <- NA
for (i in 1:nrow(bestFliers)) {
  temp <- peps[peps$Inj == bestFliers$Inj[i] &
                       peps$acetSeq == bestFliers$acetSeq[i] &
                       peps$NtermAcet_H &
                  peps$Accession2 == bestFliers$Accession2[i],]
  if (nrow(temp) == 1) {
    bestFliers$HeavyArea[i] <- temp$Area[1]
  } else if (nrow(temp) < 1) {
    bestFliers$HeavyArea[i] <- NA
  } else {
    print(nrow(temp))
    print(i)
    print(bestFliers$acetSeq[i])
  }
}


bestFliers$pctAcet <- NA
bestFliers$pctAcet[!is.na(bestFliers$Area) &
                     !is.na(bestFliers$HeavyArea)] <- 100 * 
  (bestFliers$Area[!is.na(bestFliers$Area) &
                    !is.na(bestFliers$HeavyArea)] / (bestFliers$Area[!is.na(bestFliers$Area) &
                                                                       !is.na(bestFliers$HeavyArea)] +
                                                       bestFliers$HeavyArea[!is.na(bestFliers$Area) &
                                                                         !is.na(bestFliers$HeavyArea)]))


bestFliers$pctAcet[!is.na(bestFliers$Area) &
                     is.na(bestFliers$HeavyArea)] <- 100


bestFliers$pctAcet[is.na(bestFliers$Area) &
                     !is.na(bestFliers$HeavyArea)] <- 0


CombinedTandBreps <- bestFliers %>% group_by(Accession2, gene, acetSeq) %>%
  summarise(mean_HeavyArea_WT = mean(HeavyArea[strain == "WT"], na.rm = T),
            mean_Area_WT = mean(Area[strain == "WT"], na.rm = T),
            mean_HeavyArea_DEL = mean(HeavyArea[strain == DELSTRAIN], na.rm = T),
            mean_Area_DEL = mean(Area[strain == DELSTRAIN], na.rm = T),
            mean_HeavyArea_COMP = mean(HeavyArea[strain == COMPSTRAIN], na.rm = T),
            mean_Area_COMP = mean(Area[strain == COMPSTRAIN], na.rm = T),
            mean_HeavyArea_CROSS = mean(HeavyArea[strain == "cross"], na.rm = T),
            mean_Area_CROSS = mean(Area[strain == "cross"], na.rm = T),
            count_H_WT = sum(!is.na(HeavyArea[strain == "WT"])),
            count_L_WT = sum(!is.na(Area[strain == "WT"])),
            count_H_DEL = sum(!is.na(HeavyArea[strain == DELSTRAIN])),
            count_L_DEL = sum(!is.na(Area[strain == DELSTRAIN])),
            count_H_COMP = sum(!is.na(HeavyArea[strain == COMPSTRAIN])),
            count_L_COMP = sum(!is.na(Area[strain == COMPSTRAIN])),
            count_H_CROSS = sum(!is.na(HeavyArea[strain == "cross"])),
            count_L_CROSS = sum(!is.na(Area[strain == "cross"])),
            sd_HeavyArea_WT = sd(HeavyArea[strain == "WT"], na.rm = T),
            sd_Area_WT = sd(Area[strain == "WT"], na.rm = T),
            sd_HeavyArea_DEL = sd(HeavyArea[strain == DELSTRAIN], na.rm = T),
            sd_Area_DEL = sd(Area[strain == DELSTRAIN], na.rm = T),
            sd_HeavyArea_COMP = sd(HeavyArea[strain == COMPSTRAIN], na.rm = T),
            sd_Area_COMP = sd(Area[strain == COMPSTRAIN], na.rm = T),
            sd_HeavyArea_CROSS = sd(HeavyArea[strain == "cross"], na.rm = T),
            sd_Area_CROSS = sd(Area[strain == "cross"], na.rm = T),
            mean_pctAcet_WT = mean(pctAcet[strain == "WT"], na.rm = T),
            mean_pctAcet_DEL = mean(pctAcet[strain == DELSTRAIN], na.rm = T),
            mean_pctAcet_COMP = mean(pctAcet[strain == COMPSTRAIN], na.rm = T),
            mean_pctAcet_CROSS = mean(pctAcet[strain == "cross"], na.rm = T),
            sd_pctAcet_WT = sd(pctAcet[strain == "WT"], na.rm = T),
            sd_pctAcet_DEL = sd(pctAcet[strain == DELSTRAIN], na.rm = T),
            sd_pctAcet_COMP = sd(pctAcet[strain == COMPSTRAIN], na.rm = T),
            sd_pctAcet_CROSS = sd(pctAcet[strain == "cross"], na.rm = T),
            count_pctAcet_WT = sum(!is.na(pctAcet[strain == "WT"])),
            count_pctAcet_DEL = sum(!is.na(pctAcet[strain == DELSTRAIN])),
            count_pctAcet_COMP = sum(!is.na(pctAcet[strain == COMPSTRAIN])),
            count_pctAct_CROSS = sum(!is.na(pctAcet[strain == "cross"])),
            FC_pctAcet_WTDEL = mean_pctAcet_WT / mean_pctAcet_DEL,
            FC_pctAcet_COMPDEL = mean_pctAcet_COMP / mean_pctAcet_DEL,
            pval_WTDEL = ifelse(is.character(try(t.test(pctAcet[strain == "WT"],
                                                        pctAcet[strain == DELSTRAIN])[[3]],
                                                 silent = T)),
                                NA, t.test(pctAcet[strain == "WT"],
                                           pctAcet[strain == DELSTRAIN])[[3]]),
            pval_COMPDEL = ifelse(is.character(try(t.test(pctAcet[strain == COMPSTRAIN],
                                                          pctAcet[strain == DELSTRAIN])[[3]],
                                                   silent = T)),
                                  NA, t.test(pctAcet[strain == COMPSTRAIN],
                                             pctAcet[strain == DELSTRAIN])[[3]])
    
  )


minWTandCompObs <- 2
CombinedTandBreps2 <- CombinedTandBreps[CombinedTandBreps$count_pctAcet_WT >= minWTandCompObs &
                                          CombinedTandBreps$count_pctAcet_COMP >= minWTandCompObs,] 


CombinedTandBreps2$pvalCategoryWTDEL <- CombinedTandBreps2$pval_WTDEL
CombinedTandBreps2$pvalCategoryWTDEL[CombinedTandBreps2$mean_pctAcet_WT == 100 &
                                       CombinedTandBreps2$mean_pctAcet_DEL == 0] <- 0.001

CombinedTandBreps2$RSD_pctAcet_WT <- CombinedTandBreps2$sd_pctAcet_WT / CombinedTandBreps2$mean_pctAcet_WT
for (i in 1:nrow(CombinedTandBreps2)) {
  if ((CombinedTandBreps2$RSD_pctAcet_WT[i] == 0 &
      !is.na(CombinedTandBreps2$RSD_pctAcet_WT[i]))|
      (CombinedTandBreps2$mean_pctAcet_WT[i] == 0 &
       !is.na(CombinedTandBreps2$mean_pctAcet_WT[i]))) {
    if (CombinedTandBreps2$mean_pctAcet_WT[i] == 100) {
      CombinedTandBreps2$RSD_pctAcet_WT[i] <- CombinedTandBreps2$sd_Area_WT[i] /
        CombinedTandBreps2$mean_Area_WT[i]
    } else if (CombinedTandBreps2$mean_pctAcet_WT[i] == 0) {
      CombinedTandBreps2$RSD_pctAcet_WT[i] <- CombinedTandBreps2$sd_HeavyArea_WT[i] /
        CombinedTandBreps2$mean_HeavyArea_WT[i]
    }
  }
}

CombinedTandBreps2$RSD_pctAcet_DEL <- CombinedTandBreps2$sd_pctAcet_DEL / CombinedTandBreps2$mean_pctAcet_DEL
for (i in 1:nrow(CombinedTandBreps2)) {
  if ((CombinedTandBreps2$RSD_pctAcet_DEL[i] == 0 &
      !is.na(CombinedTandBreps2$RSD_pctAcet_DEL[i])) |
      (CombinedTandBreps2$mean_pctAcet_DEL[i] == 0 &
       !is.na(CombinedTandBreps2$mean_pctAcet_DEL[i]))) {
    if (CombinedTandBreps2$mean_pctAcet_DEL[i] == 100) {
      CombinedTandBreps2$RSD_pctAcet_DEL[i] <- CombinedTandBreps2$sd_Area_DEL[i] /
        CombinedTandBreps2$mean_Area_DEL[i]
    } else if (CombinedTandBreps2$mean_pctAcet_DEL[i] == 0) {
      CombinedTandBreps2$RSD_pctAcet_DEL[i] <- CombinedTandBreps2$sd_HeavyArea_DEL[i] /
        CombinedTandBreps2$mean_HeavyArea_DEL[i]
    }
  }
}

CombinedTandBreps2$RSD_pctAcet_COMP <- CombinedTandBreps2$sd_pctAcet_COMP / CombinedTandBreps2$mean_pctAcet_COMP
for (i in 1:nrow(CombinedTandBreps2)) {
  if ((CombinedTandBreps2$RSD_pctAcet_COMP[i] == 0 &
       !is.na(CombinedTandBreps2$RSD_pctAcet_COMP[i])) |
      (CombinedTandBreps2$mean_pctAcet_COMP[i] == 0 &
       !is.na(CombinedTandBreps2$mean_pctAcet_COMP[i]))) {
    if (CombinedTandBreps2$mean_pctAcet_COMP[i] == 100) {
      CombinedTandBreps2$RSD_pctAcet_COMP[i] <- CombinedTandBreps2$sd_Area_COMP[i] /
        CombinedTandBreps2$mean_Area_COMP[i]
    } else if (CombinedTandBreps2$mean_pctAcet_COMP[i] == 0) {
      CombinedTandBreps2$RSD_pctAcet_COMP[i] <- CombinedTandBreps2$sd_HeavyArea_COMP[i] /
        CombinedTandBreps2$mean_HeavyArea_COMP[i]
    }
  }
}



CombinedTandBreps2$RSD_pctAcet_CROSS <- CombinedTandBreps2$sd_pctAcet_CROSS / CombinedTandBreps2$mean_pctAcet_CROSS
for (i in 1:nrow(CombinedTandBreps2)) {
  if ((CombinedTandBreps2$RSD_pctAcet_CROSS[i] == 0 &
       !is.na(CombinedTandBreps2$RSD_pctAcet_CROSS[i])) |
      (CombinedTandBreps2$mean_pctAcet_CROSS[i] == 0 &
       !is.na(CombinedTandBreps2$mean_pctAcet_CROSS[i]))) {
    if (CombinedTandBreps2$mean_pctAcet_CROSS[i] == 100) {
      CombinedTandBreps2$RSD_pctAcet_CROSS[i] <- CombinedTandBreps2$sd_Area_CROSS[i] /
        CombinedTandBreps2$mean_Area_CROSS[i]
    } else if (CombinedTandBreps2$mean_pctAcet_CROSS[i] == 0) {
      CombinedTandBreps2$RSD_pctAcet_CROSS[i] <- CombinedTandBreps2$sd_HeavyArea_CROSS[i] /
        CombinedTandBreps2$mean_HeavyArea_CROSS[i]
    }
  }
}



plot <- CombinedTandBreps2[CombinedTandBreps2$FC_pctAcet_WTDEL > 1 &
                             !is.na(CombinedTandBreps2$FC_pctAcet_WTDEL) &
                             CombinedTandBreps2$FC_pctAcet_COMPDEL > 1 &
                             !is.na(CombinedTandBreps2$FC_pctAcet_COMPDEL),]



plotLong <- plot %>% pivot_longer(cols = (starts_with("mean_pctAcet_")|
                                                starts_with("count_pctAcet_") |
                                                starts_with("sd_pctAcet_") |
                                            starts_with("RSD_pctAcet_")),
                                      names_to = c(".value", "strain"),
                                      names_pattern = "(.*)_pctAcet_(.*)")

plotLong[plotLong == "DEL"] <- "Δ<i>mbtK</i>"
plotLong[plotLong == "COMP"] <- "Δ<i>mbtK</i>/p<i>mbtK</i>"
plotLong[plotLong == "CROSS"] <- "Δ<i>mbtK</i>/p<i>eis</i>"

plotLong$strain <- factor(plotLong$strain, levels = c("WT", "Δ<i>mbtK</i>", "Δ<i>mbtK</i>/p<i>mbtK</i>", "Δ<i>mbtK</i>/p<i>eis</i>"))


labelDF <- read_xlsx("MMAR gene conversions for Acetylation Labels.xlsx")

for(i in 1:nrow(labelDF)) {
  accession <- labelDF$MMAR[i]
  label <- labelDF$protein_label[i]
  plotLong$gene[plotLong$Accession2 == accession] <- label
}



ggplot(plotLong) +
  geom_line(aes(x = strain, y = mean,
                group = acetSeq, color = gene),
            lwd = 1.2) +
  geom_point(aes(x = strain, y = mean,
                 group = acetSeq, fill = gene),
             stroke = 2, shape = 21, color = "black", size = 4) +
  theme_bw(base_size = 20) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_markdown()) +
  labs(y = "Percent N-terminal Acetylation",
       x = element_blank()) +
  geom_label_repel(data = plotLong[plotLong$strain == "WT",],
                   aes(x = strain, y = mean,
                       label = gene, color = gene),
                   nudge_x = -1.3, alpha = 0.9,
                   size = 5,
                   max.overlaps = 50) +
  coord_cartesian(ylim = c(0,120))

ggsave(paste0("proteins5_pub/", STRAIN, "AcetylationPlot.png"),
       width = 12, height = 14)

write.table(plot, paste0("proteins5_pub/", STRAIN, "_pct_Acet.tsv"),
            sep = "\t", row.names = F, quote = F)
