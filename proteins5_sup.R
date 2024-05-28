#set working directory####
setwd("~/04_Champion_Lab/02_N-terminal_Acetylation/3740_3692_DDABUP/DDH 3740 All LFQ/")

#import packages
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(ggrepel)
library(plotly)
library(QFeatures)
library(corrplot)
library(limma)
library(tibble)
library(imputeLCMD)

#####
#set sample type. "P" for pellet, "S" for supernatant.
sampleType <- "S"

#set pvalue cutoff. This will be B-H corrected p-values
pvalcutoff <- 0.05

#import Peaks Protein file
proteins <- read.csv("db.proteins.csv")
#remove contaminants
proteins <- proteins[!grepl("CONTAM", proteins$Accession),]
#parse out Accessions and gene names
proteins$Accession2 <- str_remove(str_extract(proteins$Accession, "^[^\\|]*\\|"), "\\|")
proteins$gene <- str_remove(str_remove(str_extract(proteins$Accession, "^[^\\|]*\\|[^\\|]*\\|"),
                                           "^[^\\|]*\\|"), "\\|")
#import to Qfeqtures, selecting area columns as quantitative, and the accessions as identifiers (row names)
Qprot <- readQFeatures(table = proteins,
                       ecol = which(grepl(paste0("Area.x[0-9]+L_", sampleType), names(proteins))),
                       fnames = "Accession2",
                       name = "raw_proteins")
Qprot #print
colnames(Qprot[['raw_proteins']]) #double check column names
#format column names
colnames(Qprot[["raw_proteins"]]) <- str_remove(colnames(Qprot[['raw_proteins']]), "Area.x[0-9]+L_")
#add metadata relating to strain, br, tr, and sample type, along with inj which has all.
Qprot$strain <- str_remove(str_extract(colnames(Qprot[["raw_proteins"]]), "\\_[^_]+$"), "_")
Qprot$bioRep <- str_remove(str_extract(colnames(Qprot[["raw_proteins"]]), "[PS][0-9]"), "[PS]")
Qprot$techRep <- str_remove_all(str_extract(colnames(Qprot[["raw_proteins"]]), "_[0-9]_"), "_")
Qprot$sampleType <- str_extract(colnames(Qprot[["raw_proteins"]]), "^[PS]")
Qprot$Inj <- colnames(Qprot[["raw_proteins"]])

colData(Qprot) #print to check

#isolate out sample type (P or S) that was set at the top
Qprot <- Qprot[Qprot$sampleType == sampleType]


#keep only Accession and gene info for each row
rowVars <- c("Accession2", "gene")
Qprot <- selectRowData(Qprot, rowVars)
unlist(rowDataNames(Qprot))

#log transform area measurements. In log2 so LFC is in log2 later.
Qprot <- addAssay(Qprot,
                  logTransform((Qprot[["raw_proteins"]]),
                               base = 2),
                  name = "log_proteins")

#normalize log values using center medians
Qprot <- addAssay(Qprot,
                  normalize(Qprot[["log_proteins"]], method = "center.median"),
                  name = "norm_proteins")

#save figure showing normalization
png("proteins5/QfeatNormAll_S.png", units = "px", width = 1200, height = 800)
par(mfrow = c(1, 3))
limma::plotDensities(assay(Qprot[[1]]), legend = F)
limma::plotDensities(assay(Qprot[[2]]), legend = F)
limma::plotDensities(assay(Qprot[[3]]), legend = F)
dev.off()

#add the same column data to the normalized protein experiment
colData(Qprot[["norm_proteins"]]) <- colData(Qprot)

#create a PCA plot
prot_pca <- Qprot[["norm_proteins"]] %>%
  filterNA() %>%
  assay() %>%
  t() %>%
  prcomp(scale = TRUE, center = TRUE)
pca_df <- as.data.frame(prot_pca$x)
pca_df$condition <- Qprot$strain
pca_df$br <- Qprot$bioRep


## Generate a PCA plot using PC1 and PC2, colored by strain, shape by br
pca_df %>%
  ggplot(aes(x = PC1, y = PC2, colour = condition, shape = br)) +
  geom_point(size = 4) +
  scale_color_brewer(palette = "Set2") +
  labs(colour = "Strain",shape = "bioRep") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(x = "PC1 (42.5 %)", y = "PC2 (26.8 %)") +
  ggtitle("Protein-level PCA plot")+
  coord_fixed(ratio = 1) +
  theme_bw()
#save
ggsave("proteins5/PCA plot for sup.png", width = 8, height = 8)

#pull out normalized proteins as a dataframe to plot correlation
prot_df <- Qprot[["norm_proteins"]] %>% assay() %>% as.data.frame()
#create correlation matrix
corr_matrixP <- cor(prot_df[, grep("^S", names(prot_df))],
                   method = "pearson",
                   use = "pairwise.complete.obs")
layout(matrix(c(1)), 1)

#save correlation plot
png("proteins5/corrplotSUP.png", width = 1000, height = 1000)
corrplot(corr_matrixP,
         method = "color",
         tl.col = "black",
         number.font = 1,
         tl.cex = 0.65,
         type = "upper")
dev.off()

#dif expression BR collapsed####
#pull out normalized proteins to collapse bioreps from technical reps, so corred DOF from BR are used
brSummary <- Qprot[["norm_proteins"]] %>% assay() %>% data.frame() %>% rownames_to_column("Accession2")

#pivot to long format pulling out BR, TR, and strain info
brSummaryL <- brSummary %>% pivot_longer(cols = which(!grepl("Accession2", names(brSummary))),
                                         names_to = c("BR", "TR", "strain"),
                                         names_pattern = "[PS](\\d+)_(\\d+)_(.*)",
                                         values_to = "Value")

#group technical reps, taking the average while ignoring NA values
brSummaryComb <- brSummaryL %>% group_by(BR, strain, Accession2) %>%
  summarise(Area = mean(Value, na.rm = T))

#pivot back to wide format taking the average Area values as the column values
brSummaryWide <- brSummaryComb %>% pivot_wider(names_from = c(BR, strain),
                                               values_from = Area,
                                               names_sep = "_")
#add sample type back to names of columns
names(brSummaryWide)[2:ncol(brSummaryWide)] <- paste0(sampleType, names(brSummaryWide)[2:ncol(brSummaryWide)])

#add back to a Qfeatures object
QprotBR <- readQFeatures(brSummaryWide,
                         ecol = 2:ncol(brSummaryWide),
                         fnames = "Accession2",
                         name = "norm_proteins_comb")

#print column names to check
colnames(QprotBR[["norm_proteins_comb"]])

#add same meta data as above
QprotBR$strain <- str_remove(str_extract(colnames(QprotBR[["norm_proteins_comb"]]), "\\_[^_]+$"), "_")
QprotBR$bioRep <- str_remove(str_extract(colnames(QprotBR[["norm_proteins_comb"]]), "[PS][0-9]"), "[PS]")
QprotBR$sampleType <- str_extract(colnames(QprotBR[["norm_proteins_comb"]]), "^[PS]")

#pull back out as single object, not necessary but makes syntax easier to read
proteins <- QprotBR[["norm_proteins_comb"]]

#set the strain as the factor of interest, and set WT as the reference
proteins$strain <- factor(QprotBR$strain)
proteins$strain <- relevel(proteins$strain, ref = "WT")

#create a linear matrix model using model.matrix to assign which columns go with which coefficients
model_design <- model.matrix(~proteins$strain)

#fit proteins to linear model with lmFit and the matrix we just created
fitted_lm <- proteins %>% 
  assay() %>%
  lmFit(design = model_design) %>% #using limma 
  # treat(lfc = 1) %>% #this would filter to things that were significantly greater or less than lfc of 1, instead of different than 0
  eBayes() # empirical bayesean technique from limma to caluclate stats using moderation of SE toward global. increases power

#makes the result easy to read in a table, sorted by most significant genes
#select each coefficient to see it related to WT.
#adjust by benjamini hochburg for adj.P.values for multiple hypothesis testing
treat_results_DEL <- topTable(fit = fitted_lm,
                          adjust.method = "BH", #BH correction
                          number = Inf, #show everything
                          coef = "proteins$straind3692", #select which comparison to make
                          sort.by = "p", confint = T) #sort by p-value (doesn't actually matter) and include 95% CIs
treat_results_DEL$Accession2 <- rownames(treat_results_DEL) #add row names based on accession


#same as above but for complement strain
treat_results_COMP <- topTable(fit = fitted_lm,
                          adjust.method = "BH",
                          number = Inf,
                          coef = "proteins$strainc3692",
                          sort.by = "p", confint = T) 
treat_results_COMP$Accession2 <- rownames(treat_results_COMP)

#same as above but cross complement strain
treat_results_CROSS <- topTable(fit = fitted_lm,
                               adjust.method = "BH",
                               number = Inf,
                               coef = "proteins$straincross",
                               sort.by = "p", confint = T) 
treat_results_CROSS$Accession2 <- rownames(treat_results_CROSS)


#add in missing and infinite values, along with genes, other info etc####
#same import as above becuase I've somewhat stupidly reused a variable name...
rawImport <- read.csv("db.proteins.csv")
rawImport <- rawImport[!grepl("CONTAM", rawImport$Accession),]
rawImport$Accession2 <- str_remove(str_extract(rawImport$Accession, "^[^\\|]*\\|"), "\\|")
rawImport$gene <- str_remove(str_remove(str_extract(rawImport$Accession, "^[^\\|]*\\|[^\\|]*\\|"),
                                       "^[^\\|]*\\|"), "\\|")

#create a dictionary linking accessions and genes
geneDict <- as.list(rawImport$gene)
names(geneDict) <- rawImport$Accession2

#pull out proteins as DF again
proteinsDF <- proteins %>% assay() %>% data.frame()

#function to add infinite FC and WT expression to df
add_infintes_and_S1_Expression_to_DF <- function(result_DF, strain1 = "WT", strain2) {

  proteinsDFs1s2 <- proteinsDF[,grepl(strain1, names(proteinsDF)) |
                                 grepl(strain2, names(proteinsDF))] #use only strains of interest
  
  proteinsDFs1s2$Accession2 <- rownames(proteinsDFs1s2) #row names as accession
  proteinsDFlong <- proteinsDFs1s2 %>% pivot_longer(1:6) #long format
  proteinsDFlong$strain <- str_remove(str_extract(proteinsDFlong$name, "_.*$"), "_") #parse strain
  proteinsDFstrains <- proteinsDFlong %>% group_by(strain, Accession2) %>%
    summarise(mean = mean(value, na.rm = T)) #group by strain and protein for mean value
  
  s1values <- proteinsDFstrains[proteinsDFstrains$strain == strain1,] #only looking at WT (strain 1 here)
  s1values$mean[is.na(s1values$mean)] <- min(s1values$mean, na.rm = T) #set NA values to minimum value
  #create dictionary linking means and accessions, for x axis of plot (WT expression levels)
  s1Dict <- as.list(s1values$mean)
  names(s1Dict) <- s1values$Accession2
  
  #pull out values containing NAs
  NAproteinsDF <- proteinsDFs1s2[rowSums(is.na(proteinsDFs1s2)) > 0 ,]
  NAproteinsLong <- NAproteinsDF %>% pivot_longer(1:6) #pivot longer 
  NAproteinsLong$strain <- str_remove(str_extract(NAproteinsLong$name, "_.*$"), "_") #parse strain
  NAproteinsStrain <- NAproteinsLong %>% group_by(Accession2) %>% #group by accessions
    summarise(means1 = mean(value[strain == strain1], na.rm = T), #pull out means for each strain
              means2 = mean(value[strain == strain2], na.rm = T),
              countNA1 = sum(is.na(value[strain == strain1])),
              countNA2 = sum(is.na(value[strain == strain2]))) #pull out NA counts for each strain
  NAproteinsStrain <- NAproteinsStrain[NAproteinsStrain$countNA1 + NAproteinsStrain$countNA2 < 6,]#remove rows without all Nas
  NAproteinsStrain$logFC <- ifelse(is.na(NAproteinsStrain$means1) |
                                     is.na(NAproteinsStrain$means2), NA,
                                   NAproteinsStrain$means2 - NAproteinsStrain$means1) #set logFC for proteins where there are values for both WT and strain
  NAproteinsStrain$logFC[is.na(NAproteinsStrain$means1)] <- Inf #set infinite for missing WT
  NAproteinsStrain$logFC[is.na(NAproteinsStrain$means2)] <- -Inf #set -infinite for missing strain
  
  #pull out only columns of interest (not actually necessary I guess...)
  NAproteinsStrain <- NAproteinsStrain[,names(NAproteinsStrain) %in% c("Accession2", "logFC")]
  NAdictlfc <- as.list(NAproteinsStrain$logFC) #create dictionary linking logfc of NA containin prots with Accessions
  names(NAdictlfc) <- NAproteinsStrain$Accession2

  #couldn't get this to work without a loop, I'm probably being silly
  #set values in the limma result DF to the lfc if they had a valid value. 
  #limma throws out rows where any value is NA. we don't do stats on these, but we want to show them on the plt.
  for(i in 1:nrow(result_DF)) {
    if(is.na(result_DF$logFC[i])) {
      if (length(NAdictlfc[[result_DF$Accession2[i]]]) > 0) {
        result_DF$logFC[i] <- NAdictlfc[[result_DF$Accession2[i]]]
      }
    }
  }
  
  
  #add gene info to result
  result_DF$gene <- sapply(result_DF$Accession2, function(Accession2) geneDict[[Accession2]])
  #add WT expression info to result (set as x axis)
  result_DF$xaxis = sapply(result_DF$Accession2, function(Accession2) s1Dict[[Accession2]])
  #set y axis variable as the logFC
  result_DF$yaxis = result_DF$logFC
  
  return(result_DF)
}

#run previous function on each result dataframe
treat_results_DEL <- add_infintes_and_S1_Expression_to_DF(treat_results_DEL,
                                                          strain1 = "WT",
                                                          strain2 = "d3692")
treat_results_COMP <- add_infintes_and_S1_Expression_to_DF(treat_results_COMP,
                                                          strain1 = "WT",
                                                          strain2 = "c3692")

treat_results_CROSS <- add_infintes_and_S1_Expression_to_DF(treat_results_CROSS,
                                                           strain1 = "WT",
                                                           strain2 = "cross")

#pull out proteins that are significant from WT vs Deletion. These will be highlighted
DELproteins <- na.omit(unique(treat_results_DEL$Accession2[treat_results_DEL$adj.P.Val < pvalcutoff]))

#Individual Plots if desired ####
colors <- c("green4", "hotpink2")
names(colors) <- c(TRUE, FALSE)
# ggplot() +
#   geom_point(data = treat_results_DEL,
#              aes(x = xaxis, y = yaxis), color = "grey50") +
#   theme_bw() +
#   geom_hline(yintercept = 0, linetype = 1) +
#   geom_hline(yintercept = c(-1,1), linetype = 2) +
#   geom_segment(data = treat_results_DEL[treat_results_DEL$adj.P.Val < pvalcutoff,],
#                aes(x = xaxis, y = CI.L, yend = CI.R), lwd = 1) +
#   geom_point(data = treat_results_DEL[treat_results_DEL$adj.P.Val < pvalcutoff,],
#              aes(x = xaxis, y = yaxis, fill = adj.P.Val < pvalcutoff), shape = 21,
#              size = 3, stroke = 2) +
#   labs(y = "Log2 Fold Change", x = "WT Relative Expression log (0 = median)",
#        title = "WT vs DEL") +
#   geom_label_repel(data = treat_results_DEL[treat_results_DEL$adj.P.Val < pvalcutoff,],
#                    aes(x = xaxis, y = yaxis, label = gene),
#                    nudge_x = 1.5, nudge_y = 4, alpha = 0.75) +
#   theme(legend.position = "none") +
#   scale_fill_manual(values = colors)
# 
# 
# DELproteins <- na.omit(unique(treat_results_DEL$Accession2[treat_results_DEL$adj.P.Val < pvalcutoff]))
# 
# 
# 
# ggplot() +
#   geom_point(data = treat_results_COMP,
#              aes(x = xaxis, y = yaxis), color = "grey50") +
#   geom_point(data = treat_results_COMP[treat_results_COMP$adj.P.Val < pvalcutoff,],
#              aes(x = xaxis, y = yaxis), color = "green2") +
#   geom_hline(yintercept = 0, linetype = 1) +
#   geom_hline(yintercept = c(-1,1), linetype = 2) +
#   geom_segment(data = treat_results_COMP[treat_results_COMP$Accession2 %in% DELproteins,],
#                aes(x = xaxis, y = CI.L, yend = CI.R), lwd = 1) +
#   geom_point(data = treat_results_COMP[treat_results_COMP$Accession2 %in% DELproteins,],
#              aes(x = xaxis, y = yaxis, fill = adj.P.Val < pvalcutoff), shape = 21,
#              size = 3, stroke = 2) +
#   theme_bw() +
#   labs(y = "Log2 Fold Change", x = "WT Relative Expression log (0 = median)",
#        title = "WT vs COMP") +
#   geom_label_repel(data = treat_results_COMP[treat_results_COMP$Accession2 %in% DELproteins,],
#                    aes(x = xaxis, y = yaxis, label = gene),
#                    nudge_x = 2, nudge_y = 2, alpha = 0.8) +
#   geom_label_repel(data = treat_results_COMP[treat_results_COMP$adj.P.Val <pvalcutoff,],
#                    aes(x = xaxis, y = yaxis, label = gene),
#                    nudge_x = 2, nudge_y = 2, color = "green4", alpha = 0.8) +
#   theme(legend.position = "none") +
#   scale_fill_manual(values = colors)


#Combined Plotting ####

#add a comparison column to differentiate WT vs what
treat_results_DEL$comparison <- "d3692"
treat_results_COMP$comparison <- "c3692"
treat_results_CROSS$comparison <-"cross"

#combine each comparison into one dataframe
combined_results <- bind_rows(treat_results_DEL, treat_results_COMP, treat_results_CROSS)
#order the comparison categories
combined_results$comparison <- factor(combined_results$comparison, levels = c("d3692",
                                                                          "c3692",
                                                                          "cross"))

#calculate corrected z score based on B-H adujsted p-value. Altman, https://doi.org/10.1136/bmj.d2090
combined_results$corrected_z <- sqrt(0.743 - (2.404*log(combined_results$adj.P.Val))) - 0.862
#calculate the one directional 95% confidence interval
combined_results$error <-  abs(combined_results$logFC /combined_results$corrected_z) * 1.96

#set the LFC limits of the plot (some 95% CI will exceed this limit)
plotLimitLFC <- 15

#plot
ggplot() +
  facet_grid(comparison~.) + #plots for each different comparison
  geom_point(data = combined_results,
             aes(x = xaxis, y = yaxis), color = "grey50", alpha = 0.8, size = 0.8) + #plot all points grey
  geom_point(data = combined_results[combined_results$adj.P.Val < pvalcutoff &
                                       !is.na(combined_results$adj.P.Val),],
             aes(x = xaxis, y = yaxis), color = "green2", alpha = 0.8, size = 0.8) + #plot significan points green
  geom_hline(yintercept = 0, linetype = 1) + #add line at LFC = 0
  geom_hline(yintercept = c(-1,1), linetype = 2) + #add lines at LFC = -1, 1
  geom_segment(data = combined_results[combined_results$Accession2 %in% DELproteins,],
               aes(x = xaxis, y = yaxis - error, yend = yaxis + error), lwd = 0.9,) + #add 95% CI lines
  geom_segment(data = combined_results[combined_results$Accession2 %in% DELproteins,],
               aes(x = xaxis, y = yaxis - error, yend = yaxis + error,
                   color = adj.P.Val < pvalcutoff), lwd = 0.4) + #there are two of these to create the color w/ black outline effect
  # geom_segment(data = combined_results[combined_results$Accession2 %in% DELproteins,],
  #              aes(x = xaxis, y = CI.L, yend = CI.R), lwd = 2, color = "red") +  #THIS WAS REMOVED, but plots the uncorrected 95% CI if that is wanted
  geom_point(data = combined_results[combined_results$Accession2 %in% DELproteins,],
             aes(x = xaxis, y = yaxis, fill = adj.P.Val < pvalcutoff), shape = 21,
             size = 3, stroke = 2) + #add points for proteins of interest
  theme_bw(base_size = 25) +
  labs(y = "Log Fold Change", x = "Normalized WT Expression log (0 = median)") +
  geom_label_repel(data = combined_results[combined_results$Accession2 %in% DELproteins,],
                   aes(x = xaxis, y = yaxis, label = gene),
                   nudge_x = 2, nudge_y = 4, alpha = 0.8, max.overlaps = 10) + #add labels to proteins of interest
  geom_label_repel(data = combined_results[combined_results$adj.P.Val <pvalcutoff &
                                             !is.na(combined_results$adj.P.Val) &
                                             !combined_results$Accession2 %in% DELproteins,],
                   aes(x = xaxis, y = yaxis, label = gene),
                   nudge_x = 2, nudge_y = 4, color = "green4", alpha = 0.8,
                   max.overlaps = 20) + #add labels to other significant proteins not on the of interest list
  theme(legend.position = "none") +
  geom_point(data = combined_results[!is.finite(combined_results$yaxis),],
              aes(x = xaxis,y = yaxis), shape = 3, size = 4, alpha = 0.8, color = "grey50") + #add infinite FC, which end up being on the axes
  scale_fill_manual(values = colors) + #set colors
  scale_color_manual(values = colors) + #set colors
  coord_cartesian(ylim = c(-plotLimitLFC, plotLimitLFC)) #set y-axis limits
  
#save
ggsave("proteins5/CombinedBayeseanBHcorrectedWTexpressionVSL2F_sup.png", width = 13, height = 18)



#save dataframe as tsv
combined_results$WTexpression <- combined_results$xaxis
write.table(combined_results, "proteins5/AnalyzedResults_SUP.tsv",
            sep = "\t", quote = F, row.names = F)



#pull out mbt and ecc proteins in two lists. Variable name is not quite re
mtbs <- na.omit(unique(combined_results$gene[grepl("^mbt", combined_results$gene)]))
eccs <- na.omit(unique(combined_results$gene[grepl("^ecc", combined_results$gene)]))



#set protein of interest to desired list, comment out the other one.
DELproteins <- mtbs
DELproteins <- eccs


ggplot() + # plot, same as above.
  facet_grid(comparison~.) +
  geom_point(data = combined_results[!is.finite(combined_results$yaxis),],
             aes(x = xaxis,y = yaxis), shape = 3, size = 4, alpha = 0.8, color = "grey50") +
  geom_point(data = combined_results,
             aes(x = xaxis, y = yaxis), color = "grey50", alpha = 0.8, size = 0.8) +
  geom_point(data = combined_results[combined_results$adj.P.Val < pvalcutoff &
                                       !is.na(combined_results$adj.P.Val),],
             aes(x = xaxis, y = yaxis), color = "green2", alpha = 0.8, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 1) +
  geom_hline(yintercept = c(-1,1), linetype = 2) +
  geom_segment(data = combined_results[combined_results$gene %in% DELproteins,],
               aes(x = xaxis, y = yaxis - error, yend = yaxis + error), lwd = 0.9,) +
  geom_segment(data = combined_results[combined_results$gene %in% DELproteins,],
               aes(x = xaxis, y = yaxis - error, yend = yaxis + error,
                   color = adj.P.Val < pvalcutoff), lwd = 0.4)+
  # geom_segment(data = combined_results[combined_results$Accession2 %in% DELproteins,],
  #              aes(x = xaxis, y = CI.L, yend = CI.R), lwd = 2, color = "red") +
  geom_point(data = combined_results[combined_results$gene %in% DELproteins,],
             aes(x = xaxis, y = yaxis, fill = adj.P.Val < pvalcutoff), shape = 21,
             size = 3, stroke = 2) +
  theme_bw(base_size = 25) +
  labs(y = "Log Fold Change", x = "Normalized WT Expression log (0 = median)") +
  geom_label_repel(data = combined_results[combined_results$gene %in% DELproteins,],
                   aes(x = xaxis, y = yaxis, label = gene),
                   nudge_x = 0.5, nudge_y = 4, alpha = 0.8, max.overlaps = 20) +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  coord_cartesian(ylim = c(-plotLimitLFC, plotLimitLFC))

#save, use the line that corresponds to correct protein of interest list (mbt or ecc)
ggsave("proteins5/CombinedBayeseanBHcorrectedWTexpressionVSL2F_eccs_SUP.png", width = 13, height = 18)
#ggsave("proteins5/CombinedBayeseanBHcorrectedWTexpressionVSL2F_mbts_SUP.png", width = 13, height = 18)


#esxs
esxs <- c("MMAR_5450", "MMAR_5457", "MMAR_5448", "MMAR_5449", "MMAR_4166", "MMAR_2894", "MMAR_5440", "MMAR_5439", "MMAR_5453", "MMAR_5455")
DELproteins <- esxs


ggplot() + # plot, same as above.
  facet_grid(comparison~.) +
  geom_point(data = combined_results[!is.finite(combined_results$yaxis),],
             aes(x = xaxis,y = yaxis), shape = 3, size = 4, alpha = 0.8, color = "grey50") +
  geom_point(data = combined_results,
             aes(x = xaxis, y = yaxis), color = "grey50", alpha = 0.8, size = 0.8) +
  geom_point(data = combined_results[combined_results$adj.P.Val < pvalcutoff &
                                       !is.na(combined_results$adj.P.Val),],
             aes(x = xaxis, y = yaxis), color = "green2", alpha = 0.8, size = 0.8) +
  geom_hline(yintercept = 0, linetype = 1) +
  geom_hline(yintercept = c(-1,1), linetype = 2) +
  geom_segment(data = combined_results[combined_results$Accession2 %in% DELproteins,],
               aes(x = xaxis, y = yaxis - error, yend = yaxis + error), lwd = 0.9,) +
  geom_segment(data = combined_results[combined_results$Accession2 %in% DELproteins,],
               aes(x = xaxis, y = yaxis - error, yend = yaxis + error,
                   color = adj.P.Val < pvalcutoff), lwd = 0.4)+
  # geom_segment(data = combined_results[combined_results$Accession2 %in% DELproteins,],
  #              aes(x = xaxis, y = CI.L, yend = CI.R), lwd = 2, color = "red") +
  geom_point(data = combined_results[combined_results$Accession2 %in% DELproteins,],
             aes(x = xaxis, y = yaxis, fill = adj.P.Val < pvalcutoff), shape = 21,
             size = 3, stroke = 2) +
  theme_bw(base_size = 25) +
  labs(y = "Log Fold Change", x = "Normalized WT Expression log (0 = median)") +
  geom_label_repel(data = combined_results[combined_results$Accession2 %in% DELproteins,],
                   aes(x = xaxis, y = yaxis, label = gene),
                   nudge_x = 0.5, nudge_y = 4, alpha = 0.8, max.overlaps = 20) +
  theme(legend.position = "none") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  coord_cartesian(ylim = c(-plotLimitLFC, plotLimitLFC))

ggsave("proteins5/CombinedBayeseanBHcorrectedWTexpressionVSL2F_esx_SUP.png", width = 13, height = 18)


