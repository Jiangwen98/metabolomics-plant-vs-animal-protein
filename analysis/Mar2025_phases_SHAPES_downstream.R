# ---------------------------------------------------------------------------
# NOTE ON THIS PUBLIC COPY
# Trial participant identifiers have been REDACTED and replaced with
# <..._REDACTED...> placeholders. Absolute local paths have been made relative.
# No participant data is included in this repository. See README.md.
# ---------------------------------------------------------------------------

#* To-do
#* 1. add in other possible contrast dimensions - refer to the first meeting slides
#* 2. b52(? to double check ion-chro to confirm) - W0 data of SHP039 is abnormal (probably due to assay problem), need to pay attention

library(tidyverse)
library(dplyr)

### Add-in group information
#** Remove BQC first
quant_af <- quant_af[-which(quant_af$type == "BQC"),]
#** read in grouping information
grp <- read.delim("SHAPES_group_info.txt", 
                  as.is = T, check.names = F)

#** add in group information
quant_info <- NULL

for(i in 1:nrow(quant_af)){
  row <- quant_af[i,]
  fname <- row$filename
  name <- substr(fname, start = 1, stop = unlist(gregexpr("\\.", fname))-1)
  for(j in 1:nrow(grp)){
    if(grp[j,1] == name){
      row <- add_column(row, grp[j,3:5], .after = "type")
    }
  }
  rm(j)
  
  quant_info <- rbind(quant_info, row)
}
rm(i, row, fname, name)

#** remove the "type" column
quant_info <- quant_info[-3]
quant_lmm <- quant_info #LMM can handle missing data

### Data cleaning
#** remove the patients without W8 data
quant_info <- quant_info[-which(quant_info$`Patient ID` == "<SHAPES_SUBJECT_REDACTED_1>"| #PBMD
                                quant_info$`Patient ID` == "<SHAPES_SUBJECT_REDACTED_2>"| #ABMD
                                quant_info$`Patient ID` == "<SHAPES_SUBJECT_REDACTED_3>"| #PBMD
                                quant_info$`Patient ID` == "<SHAPES_SUBJECT_REDACTED_4>"| #ABMD
                                quant_info$`Patient ID` == "<SHAPES_SUBJECT_REDACTED_5>"| #PBMD
                                quant_info$`Patient ID` == "<SHAPES_SUBJECT_REDACTED_6>"),] #PBMD
#### Start statistical analysis ####
### Baseline difference comparison
quant_bl <- quant_info[which(quant_info$Week == "W0"),]
quant_bl <- quant_bl[-c(1, 2, 5)]
Res_bl <- NULL

for (i in c(3:ncol(quant_bl))) {
  
  crt_data <- quant_bl[c(1, 2, i)]
  P <- subset(crt_data, Treatment == "PBMD")
  A <- subset(crt_data, Treatment == "ABMD")
  
  ###negative data correction reference:
  ###http://younglab.wi.mit.edu/chromatin/foldchange.html
  Test <- wilcox.test(unlist(P[3]), unlist(A[3]))
  FC <- median(as.numeric(unlist(P[3]))) / median(as.numeric(unlist(A[3])))
  log2FC <- log(FC, 2)
  
  MAD_P <- mad(unlist(P[3]))
  MAD_A <- mad(unlist(A[3]))
  Med_P <- median(unlist(P[3]))
  Med_A <- median(unlist(A[3]))
  p <- Test$p.value
  
  Row <- c(colnames(quant_bl)[i], MAD_P, MAD_A, Med_P, Med_A, FC, log2FC, p)
  Res_bl <- rbind(Res_bl, Row)
  
  rm(i, crt_data, P, A, Test, FC, log2FC, MAD_P, MAD_A, Med_P, Med_A, p, Row)
}


Res_bl <- as.data.frame(Res_bl)
colnames(Res_bl) <- c("Analyte", "MAD_PBMD", "MAD_ABMD", "Median_PBMD", "Median_ABMD", "FC", "log2FC", "p.val")
Res_bl$Corrected_p <- p.adjust(unlist(Res_bl[8]), method = "BH")
Res_bl$'Neg_log10P.adj' <- -log(unlist(Res_bl[9]), 10)

write.csv(Res_bl, file = "otpt/thesis_bl_p&FC_Auto11[part00].csv", row.names = F)


### (A8-A0)/A0 vs (P8-P0)/P0 comparison
#** calculate (A8-A0)/A0 and (P8-P0)/P0
quant_patient <- quant_info %>% group_by(`Patient ID`) %>% group_split()

bl_crt_diff <- NULL # means baseline-corrected-differences
for (i in 1:length(quant_patient)){
  single_patient <- quant_patient[[i]]
  single_diff <- single_patient[1, c(3,4)]
  
  for(j in 6:ncol(single_patient)){
    this_w0 <- single_patient[[which(single_patient$Week == "W0"), j]]
    this_w8 <- single_patient[[which(single_patient$Week == "W8"), j]]
    this_diff <- (this_w8 - this_w0) / this_w0
    single_diff <- cbind(single_diff, this_diff)
    colnames(single_diff)[[ncol(single_diff)]] <- colnames(single_patient[j])
  }
  rm(j, this_w0, this_w8, this_diff)
  
  bl_crt_diff <- rbind(bl_crt_diff, single_diff)
}
rm(i, single_patient, single_diff)

write.csv(bl_crt_diff, file = "otpt/thesis_BaselineCorrected_8w_difference_Auto11[part00].csv", row.names = F)

### Subset quant_info and bl_crt_diffe for significant metaboliets after manual peak review for Xinyi to do correlation analysis ###
sig.mtbl <- read.delim("sig.mtbl.aft.peak review.txt", as.is = T, check.names = F, header = F)[[1]]
quant.sig <- quant_info[c(3:5, which(colnames(quant_info) %in% sig.mtbl))]
bl_crt.sig <- bl_crt_diff[c(1:2, which(colnames(bl_crt_diff) %in% sig.mtbl))]

# write.csv(quant.sig, file = "otpt/quant table of significant metabolites.csv", row.names = F)
# write.csv(bl_crt.sig, file = "otpt/baseline corrected difference of significant metabolites.csv", row.names = F)
### Subset quant_info and bl_crt_diffe for significant metaboliets after manual peak review for Xinyi to do correlation analysis ###




# #### pca and pls-da
# library(mixOmics)
# MR_pca_df <- MR_pls_df <- within(bl_crt_diff, Treatment <- factor(Treatment, levels = c("ABMD","PBMD"), ordered = T))
# # #pca
# pca.mr <- prcomp(as.matrix(sapply(MR_pca_df[, -c(1:2)], as.numeric)), center = TRUE, scale. = TRUE)
# pca.scores <- as.data.frame(pca.auc$x)
# pca.scores$Treatment <- MR_pca_df$Treatment
# 
# # pca.2 <- pca(as.matrix(sapply(MR_pca_df[, -c(1:2)], as.numeric)), ncomp = 2)
# # plotIndiv(pca.2)
# # plotLoadings(pca.2)
# 
# pp <- ggplot(pca.scores, aes(x = PC1, y = PC2, color = Treatment)) +
#   geom_point(size = 3) +
#   theme_minimal() +
#   labs(title = "PCA Plot",
#        x = paste0("PC1 (", round(100*summary(pca.auc)$importance[2,1], 1), "%)"),
#        y = paste0("PC2 (", round(100*summary(pca.auc)$importance[2,2], 1), "%)")) +
#   scale_color_manual(values = c("ABMD" = "red", "PBMD" = "blue"))
# pp #kind of meaningless...
# 
# ###pls-da
# pls.mr <- plsda(X = as.matrix(sapply(MR_pls_df[, -1], as.numeric)), Y = MR_pls_df$Treatment, ncomp = 2, scale = T)
# 
# plotIndiv(pls.mr, ind.names = F, ellipse = TRUE, legend = TRUE) # this is just to check the comp1 and comp2
# 
# p <- plotIndiv(pls.mr, ind.names = F, ellipse = TRUE, legend = TRUE, 
#                title = "", 
#                X.label = "Comp1 (5%)", Y.label = "Comp2 (6%)",
#                col = c("dodgerblue4", "plum"),
#                point.lwd = 0.75, cex = 1.8, style = "ggplot2")$graph + 
#   theme_classic() +
#   # ggtitle(paste0("PLS-DA of Metabolome")) + 
#   theme(plot.title = element_text(size = 17, hjust = 0.5, face = 'bold'),
#         axis.title.x = element_text(size = 12, face = 'bold'),
#         axis.title.y = element_text(size = 12, face = 'bold'),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         legend.title = element_text(size = 12, face = 'bold'),
#         legend.text = element_text(size = 12, face = 'bold'),
#         legend.position = 'right',
#         strip.background = element_blank())
# p
# 
# plotLoadings(pls.mr, contrib = "max", method = "median")
# 
# #### pca and pls-da













#** (A8-A0)/A0 vs. (P8-P0)/P0 contrast - Mann-witney U test (fold change hard to show as some in negative value)
Res <- NULL

for (i in c(3:ncol(bl_crt_diff))) {
  
  crt_data <- bl_crt_diff[c(1, 2, i)]
  P <- subset(crt_data, Treatment == "PBMD")
  A <- subset(crt_data, Treatment == "ABMD")
  
  ###negative data correction reference:
  ###http://younglab.wi.mit.edu/chromatin/foldchange.html
  MinP <- min(P[3])
  MinA <- min(A[3])
  Minz <- min(MinP, MinA)
  P_corrected <- P[3] - Minz + 5
  A_corrected <- A[3] - Minz + 5
  Test <- wilcox.test(unlist(P_corrected), unlist(A_corrected))
  FC <- median(as.numeric(unlist(P_corrected))) / median(as.numeric(unlist(A_corrected)))
  log2FC <- log(FC, 2)
 

  
  MAD_P <- mad(unlist(P[3]))
  MAD_A <- mad(unlist(A[3]))
  Med_P <- median(unlist(P[3]))
  Med_A <- median(unlist(A[3]))
  p <- Test$p.value
  
  Row <- c(colnames(bl_crt_diff)[i], MAD_P, MAD_A, Med_P, Med_A, FC, log2FC, p)
  Res <- rbind(Res, Row)
 
  rm(i, crt_data, P, A, Test, FC, log2FC, MAD_P, MAD_A, Med_P, Med_A, p, Row) 
}

Res <- as.data.frame(Res)
colnames(Res) <- c("Analyte", "MAD_PBMD", "MAD_ABMD", "Median_PBMD", "Median_ABMD", "FC", "log2FC", "p.val")
Res$Corrected_p <- p.adjust(unlist(Res[8]), method = "BH")
Res$'Neg_log10p.adj' <- -log(unlist(Res[9]), 10)

write.csv(Res, file = "otpt/thesis_bl_crt_p&FC_Auto11[part00].csv", row.names = F)


#### Using linear-mixed model to compare PBMD vs. ABMD acorrding to gpt's advice, since the mann-witney U test giving out quite less analytes with significant p.adj.
library(lme4)
library(lmerTest)
library(stats)
library(car)

## create log2-transformed data
quant_lmm[6:ncol(quant_lmm)] <- lapply(quant_lmm[6:ncol(quant_lmm)], as.numeric)
quant_lmm <- quant_lmm[-(1:2)]
quant_lmm <- quant_lmm %>% mutate(across(where(is.numeric), ~ log2(.)))

### also replace those -Inf in the outcome with NA
infinite_indices <- which(sapply(quant_lmm, is.infinite), arr.ind = TRUE)
View(infinite_indices)
# ### (1,1)
# quant_lmm[32,156] <- quant_lmm[45,156] <- NA
# ### (1,1)
# ### scaling and centering the data, to see if result has any change - no change to p-values
# quant_lmm[, 6:ncol(quant_lmm)] <- scale(quant_lmm[, 6:ncol(quant_lmm)], center = T, scale = T)

## Filling the rows contain missing W8 samples... so LMM can handle with the missing samples
#"<SHAPES_SUBJECT_REDACTED_1>","<SHAPES_SUBJECT_REDACTED_2>", "<SHAPES_SUBJECT_REDACTED_3>", "<SHAPES_SUBJECT_REDACTED_4>", "<SHAPES_SUBJECT_REDACTED_5>", "<SHAPES_SUBJECT_REDACTED_6>"
a <- c("<SHAPES_SUBJECT_REDACTED_1>", "PBMD", "W8", rep(NA, (ncol(quant_lmm)-3)))
b <- c("<SHAPES_SUBJECT_REDACTED_2>", "ABMD", "W8", rep(NA, (ncol(quant_lmm)-3)))
c <- c("<SHAPES_SUBJECT_REDACTED_3>", "PBMD", "W8", rep(NA, (ncol(quant_lmm)-3)))
d <- c("<SHAPES_SUBJECT_REDACTED_4>", "ABMD", "W8", rep(NA, (ncol(quant_lmm)-3)))
e <- c("<SHAPES_SUBJECT_REDACTED_5>", "PBMD", "W8", rep(NA, (ncol(quant_lmm)-3)))
f <- c("<SHAPES_SUBJECT_REDACTED_6>", "PBMD", "W8", rep(NA, (ncol(quant_lmm)-3)))

quant_lmm <- rbind(quant_lmm, a,b,c,d,e,f)


### Setting reference group, ABMD for treatment and W0 for Week
quant_lmm$Treatment <- relevel(factor(quant_lmm$Treatment), ref = "ABMD")
quant_lmm$Week <- relevel(factor(quant_lmm$Week), ref = "W0")
quant_lmm[4:ncol(quant_lmm)] <- lapply(quant_lmm[4:ncol(quant_lmm)], as.numeric)

colnames(quant_lmm)[[1]] <- "Patient.ID"
### Fitting LMM -  copying from my old codes
Interact.e <- mtbl.name <- Treatment.p <- Week.p <- Interact.p <- vector()
re.tab <- sum.tab <- data.frame(matrix(ncol = 6, nrow = 0))
p.tab <- data.frame(matrix(ncol = 0, nrow = (ncol(quant_lmm)-3))) ### 165 is the amount of analytes.

for (m in 4:ncol(quant_lmm)){
  ctr.name <- colnames(quant_lmm)[[m]]
  mtbl.name <- append(mtbl.name, ctr.name)
  otc.name <- paste0('`', ctr.name, '`')
  tab.col <- rep(ctr.name, 4)
  re.tab.col <- rep(ctr.name, 2)
  
  ctr.df <- quant_lmm[, c(1:3, m)]
  ctr.lm <- lmer(eval(parse(text = otc.name)) ~ Treatment + Week + Treatment:Week + (1|Patient.ID), ctr.df)
  
  
  ctr.sum <- cbind(data.frame(coef(summary(ctr.lm))), tab.col)
  sum.tab <- rbind(sum.tab, ctr.sum)
  
  ctr.re <- cbind(as.data.frame(VarCorr(ctr.lm)), re.tab.col)
  re.tab <- rbind(re.tab, ctr.re)
  
  # Treatment.p <- append(Treatment.p, data.frame(coef(summary(ctr.lm)))$`Pr...t.`[[2]]) ###!!!!!!should use anova
  # Week.p <- append(Week.p, data.frame(coef(summary(ctr.lm)))$`Pr...t.`[[3]])
  # Interact.p <- append(Interact.p, data.frame(coef(summary(ctr.lm)))$`Pr...t.`[[4]])
  
  Treatment.p <- append(Treatment.p, anova(ctr.lm)[which(rownames(anova(ctr.lm)) == "Treatment"), 6]) ###!!!!!!should use anova
  Week.p <- append(Week.p, anova(ctr.lm)[which(rownames(anova(ctr.lm)) == "Week"), 6])
  Interact.p <- append(Interact.p, anova(ctr.lm)[which(rownames(anova(ctr.lm)) == "Treatment:Week"), 6])
  
  Interact.e <- append(Interact.e, data.frame(coef(summary(ctr.lm)))$`Estimate`[[4]])
  
  print(ctr.name)
  
  rm(m, ctr.name, otc.name, tab.col, re.tab.col, ctr.df, ctr.lm, ctr.sum, ctr.re)
}
p.tab <- cbind(mtbl.name, Treatment.p, Week.p, Interact.p)

colnames(sum.tab) <- c("Estimate", "Std Error", "df", "t-value", "p-value", "biomarker")
colnames(p.tab) <- c("biomarker", "p-value_Meal", "p-value_Week", "p-value_Interaction")
colnames(re.tab)<- c("group", "var1", "var2", "variance", "Std Dev", "biomarker")

###1. [full dataset]LMM output with indicating random effect
write.csv(sum.tab, file = paste0("otpt/Thesis_ANOVA_LMM_Summary_auto11[part00]", ".csv"), row.names = T)
write.csv(p.tab, file = paste0("otpt/Thesis_ANOVA_LMM_[original]pval_auto11[part00]", ".csv"), row.names = F)
write.csv(re.tab, file = paste0("otpt/Thesis_ANOVA_LMM_Random_effect_var_auto11[part00]", ".csv"), row.names = T)

rm(mtbl.name, re.tab)

#### Modeling output p-values adjustment
library(qvalue)
if(!file.exists("otpt/pval_hist")){dir.create("otpt/pval_hist")}
if(!file.exists("otpt/pval_adj")){dir.create("otpt/pval_adj")}

p.tab <- as.data.frame(p.tab)
BH_pval_df <- qvalue_df <- p.tab[1]

for(i in 2:ncol(p.tab)){
  
  ##checking original p-value hist curve
  name <- colnames(p.tab)[[i]]
  pdf(file = paste0("otpt/pval_hist/Histogram_", name, '.pdf'), width = 11, height = 8)
  hist(x = as.numeric(p.tab[[i]]), main = paste0("Histogram of ", name), xlab = "p-values")
  dev.off()
  
  ##doing BH and q-value method for p-value adjustment
  col_suf <- str_remove(name, "^p-value_")
  BH_pval <- p.adjust(p = as.numeric(p.tab[[i]]), method = "BH") # or use qvalue(.., pi0 = 1)$qvalues
  qvalue <- qvalue(p = as.numeric(p.tab[[i]]))$qvalues
  
  BH_pval_df <- cbind(BH_pval_df, BH_pval)
  qvalue_df <- cbind(qvalue_df, qvalue)
  colnames(BH_pval_df)[[ncol(BH_pval_df)]] <- paste0("BH_pval_", col_suf)
  colnames(qvalue_df)[[ncol(qvalue_df)]] <- paste0("qvalue_", col_suf)
  
  ##graphic comparing two p-value adjustment method
  pdf(file = paste0("otpt/pval_adj/BH_vs_qvalue", "_", name, ".pdf"), width = 8, height = 8)
  plot(x = as.numeric(p.tab[[i]]), y = qvalue, xlab = "p-values", ylab = "q-values/BH adjusted p-values",
       main = paste0("BH vs q-value of ", name), cex = .1, col = 2, ylim = c(0,1), xlim = c(0,1))
  points(x = as.numeric(p.tab[[i]]), y = BH_pval, cex = .1, col = 4)
  legend("bottomright", c("q-value", "BH adjusted p-value"), pch = 19, col = c(2, 4))
  abline(a = 0, b = 1, col = 1)
  dev.off()
  
  rm(col_suf, BH_pval, qvalue, name)
}

write.csv(x = BH_pval_df, file = paste0("otpt/BH_adj_pvals_auto11[part00]", ".csv"), row.names = F)
write.csv(x = qvalue_df, file = paste0("otpt/qvalues_auto11[part00]", ".csv"), row.names = F)

##### volcano plots and box-plots of two overall comparison with adjusted pvalues added #####
### Preparing the input data for volcano plots - format: c(Name(metabolite's name), x(effect size), y(p-value))
##activate the function plot_volcano1 in function script first

vol.df <- BH_pval_df

Treatment_ES <- sum.tab[which(grepl("TreatmentPBMD$", rownames(sum.tab)) | grepl("TreatmentPBMD[0-9]+", rownames(sum.tab))), 1]
Week_ES <- sum.tab[which(grepl("^WeekW$", rownames(sum.tab)) | grepl("^WeekW[0-9]+", rownames(sum.tab))), 1]
Interaction_ES <- sum.tab[which(grepl("\\:", rownames(sum.tab))), 1]

vol.df <- add_column(vol.df, Treatment_ES, .before = "BH_pval_Meal") %>% add_column(Week_ES, .before = "BH_pval_Week") %>% add_column(Interaction_ES, .before = "BH_pval_Interaction")

write.csv(vol.df, file = paste0("otpt/LMM_BHpval_ES_auto11[part00]", ".csv"), row.names = F)

### one by one plot and adjusting...same as the one in code testing...###
Single_df <- Res[c(1, 7, 9)]##format: c(Name(metabolite's name), x(FC), y(p-value))
Single_df[2:3] <- as.numeric(unlist(Single_df[2:3]))
volc_title <- "Contrast of δMR (PBMD vs. ABMD)"
x_title <- "Log2FC"

p0 <- 
  plot_volcano1(Single_df, pval.threshold = 1, volc.title = volc_title, x.title = x_title, ####treshold change to p.adj<0.1!!!!
                y.title = '-lg(BH-adjusted p-values)', toLog2FC = F, toNegLog10Pval = T, 
                # x_limit = c(-1, 1.6),
                # y_limit = c(0, 3),
                 fc.threshold.left = 0, fc.threshold.right = 0,
                text_size = 3, repel_labels = T, dot_labels = T)

ggsave(filename = paste0("otpt/volcano of ", volc_title, ".svg")
       , plot = p0, units = "in",  width = 11, height = 8
       , device = "svg")

dev.off()

# ###temp code
# first <- unique(read.delim("p.adj.less.0.1.txt", 
#                     as.is = T, check.names = F, header = F))
# second <- unique(read.delim("p.ori.less.0.05.txt",
#                      as.is = T, check.names = F, header = F))
# to.get <- unlist(unique(read.delim("to.get.RT.txt", 
#                            as.is = T, check.names = F, header = F)))
# to.get.df <- quant_bf[1:3, which(colnames(quant_bf) %in% to.get)]
# rownames(to.get.df) <- c("Q1", "Q3", "detected_RT")
# 
# write.csv(to.get.df, file = paste0("otpt/Q1Q3RT", ".csv"), row.names = F)
# ###temp code

### Boxplot for significant metabolites of PBMD vs. ABMD

