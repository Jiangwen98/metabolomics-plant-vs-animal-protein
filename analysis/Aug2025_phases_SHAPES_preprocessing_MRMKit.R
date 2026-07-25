# ---------------------------------------------------------------------------
# NOTE ON THIS PUBLIC COPY
# Trial participant identifiers have been REDACTED and replaced with
# <..._REDACTED...> placeholders. Absolute local paths have been made relative.
# No participant data is included in this repository. See README.md.
# ---------------------------------------------------------------------------

library(tidyverse)
library(dplyr)

getwd()
setwd(".")  # path made relative for public copy
#* here I store my codes in cloud storage, so must use "~" for home directory 
#* cos the specific name of home dir may differ from work computer to my pc
if(!file.exists("otpt")){dir.create("otpt")}

###*** Read in quant table
# ### 1.(1,1)
# # quant_bf <- read.delim("mrmkit_output/2025-08-13-12-46-03 Jiangwen0513_2/quanttables/plate123_odd_even_quant_table.txt", 
# #                   as.is = T, check.names = F)
# ### 1.(1,1)
# 
# ### 2.(0,0)
# quant_bf <- read.delim("mrmkit_output/2025-08-13-13-24-50 Jiangwen0513_4/quanttables/plate123_odd_even_quant_table.txt", 
#                        as.is = T, check.names = F)
# # remove the istd columns from (0,0)
# istd_to_remove <- read.delim("istd_to_remove.txt", as.is = T, check.names =F)[[1]]
# quant_bf <- quant_bf[-which(colnames(quant_bf) %in% istd_to_remove)]
# ### 2.(0,0)

### 3.read-in (1,1), then replace those analyes' column whose istd is arginine istd or glutamic acid istd with (0,0)
quant_bf <- read.delim("mrmkit_output/2025-08-13-12-46-03 Jiangwen0513_2/quanttables/plate123_odd_even_quant_table.txt", 
                                          as.is = T, check.names = F)
mtbl_to_replace <- read.delim("mtbl_to_replace.txt", as.is = T, check.names =F)[[1]]
quant_bf <- quant_bf[-which(colnames(quant_bf) %in% mtbl_to_replace)]
quant_00 <- read.delim("mrmkit_output/2025-08-13-13-24-50 Jiangwen0513_4/quanttables/plate123_odd_even_quant_table.txt", 
                                               as.is = T, check.names = F)
to_replace_with <- quant_00[which(colnames(quant_00) %in% mtbl_to_replace)]
quant_bf <- cbind(quant_bf, to_replace_with)
### 3.read-in (1,1), then replace those analyes' column whose istd is arginine istd or glutamic acid istd with (0,0)

###*** Read in S/N from raw data
## 1.(1,1)
# sn <- read.delim("mrmkit_output/2025-08-13-10-20-08 Jiangwen0513_5/quanttables/plate123_odd_even_quant_table.txt", as.is = T, check.names = F)[10, ]
## 1.(1,1)

# ## 2.(0,0)
# quant_sn <- read.delim("mrmkit_output/2025-08-13-10-45-59 Jiangwen0513_7/quanttables/plate123_odd_even_quant_table.txt", as.is = T, check.names = F)
# quant_sn <- quant_sn[-which(colnames(quant_sn) %in% istd_to_remove)]
# sn <- quant_sn[10, ]
# ## 2.(0,0)

## 3.read-in (1,1), then replace those analyes' column whose istd is arginine istd or glutamic acid istd with (0,0)
sn_quant_11 <- read.delim("mrmkit_output/2025-08-13-10-20-08 Jiangwen0513_5/quanttables/plate123_odd_even_quant_table.txt", as.is = T, check.names = F)
sn_quant_11 <- sn_quant_11[-which(colnames(sn_quant_11) %in% mtbl_to_replace)]
sn_11 <- sn_quant_11[10, ]

sn_quant_00 <- read.delim("mrmkit_output/2025-08-13-10-45-59 Jiangwen0513_7/quanttables/plate123_odd_even_quant_table.txt", as.is = T, check.names = F)
sn_quant_00 <- sn_quant_00[which(colnames(sn_quant_00) %in% mtbl_to_replace)]
sn_00 <- sn_quant_00[10, ]
sn <- cbind(sn_11, sn_00)
## 3.read-in (1,1), then replace those analyes' column whose istd is arginine istd or glutamic acid istd with (0,0)

sn <- sn[-(2:3)]

# Now filter out using S/N, keep those with S/N >= 5
sn_pass <- colnames(sn[-1])[which(as.numeric(sn[-1]) >= 5)] ### must use colnamesof !!!sn[-1]!!! will I retrive the correct analytes!!!
### ^^^ ibatch correction won't affect S/N of specific batches, so these S/N won't change after batch corretcion.

#** QC metrics (cov and d-ratio) filtering, keep those with CoV <= 30 OR D-ratio <= 50
CoV_metric <- quant_bf[7, ][-(2:3)]
D_ratio_metric <- quant_bf[9, ][-(2:3)]

cov_pass <- colnames(CoV_metric[-1])[which(CoV_metric[-1] <= 30)]
dratio_pass <- colnames(D_ratio_metric[-1])[which(D_ratio_metric[-1] <= 50)]

analytes_pass <- intersect(union(cov_pass, dratio_pass), sn_pass)

quant_af <- quant_bf[-(1:10), ]
quant_af <- quant_af[, c(1:3, which(colnames(quant_af) %in% analytes_pass))]
