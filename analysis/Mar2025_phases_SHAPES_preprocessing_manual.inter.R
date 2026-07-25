# ---------------------------------------------------------------------------
# NOTE ON THIS PUBLIC COPY
# Trial participant identifiers have been REDACTED and replaced with
# <..._REDACTED...> placeholders. Absolute local paths have been made relative.
# No participant data is included in this repository. See README.md.
# ---------------------------------------------------------------------------

##### To research and to read the following:
### 1. cov, d-ratio after batch correction?
### 2. use mad instead of cov and d-ratio?
### 4. Try different batch correction methods: median average, grand median, ...

library(tidyverse)
library(dplyr)

getwd()
setwd(".")  # path made relative for public copy
#* here I store my codes in cloud storage, so must use "~" for home directory 
#* cos the specific name of home dir may differ from work computer to my pc
if(!file.exists("otpt")){dir.create("otpt")}

### Read in 6 quant tables (let's start with the (1,1) one)
p1e <- read.delim("mrmkit_output/2025-08-06-18-19-59 Jiangwen0513_2/quanttables/Plate1_EVEN_quant_table.txt", 
                  as.is = T, check.names = F)
p1e$batch <- 1
### Remove missing compounds that missing in Plate 1 ODD batch...
p1e <- p1e[-c(which(colnames(p1e) %in% c("3_aminoisobutanoate_a", "3_aminoisobutanoate_b", "alanine_b")))]
### Remove missing compounds that missing in Plate 1 ODD batch...


p1o <- read.delim("mrmkit_output/2025-08-06-18-19-59 Jiangwen0513_2/quanttables/Plate1_ODD_quant_table.txt", 
                  as.is = T, check.names = F)
p1o$batch <- 2


p2e <- read.delim("mrmkit_output/2025-08-06-18-19-59 Jiangwen0513_2/quanttables/Plate2_EVEN_quant_table.txt", 
                  as.is = T, check.names = F)
p2e$batch <- 3
### Remove missing compounds that missing in Plate 1 ODD batch...
p2e <- p2e[-c(which(colnames(p2e) %in% c("3_aminoisobutanoate_a", "3_aminoisobutanoate_b", "alanine_b")))]
### Remove missing compounds that missing in Plate 1 ODD batch...


p2o <- read.delim("mrmkit_output/2025-08-06-18-19-59 Jiangwen0513_2/quanttables/Plate2_ODD_quant_table.txt", 
                  as.is = T, check.names = F)
p2o$batch <- 4
### Remove missing compounds that missing in Plate 1 ODD batch...
p2o <- p2o[-c(which(colnames(p2o) %in% c("3_aminoisobutanoate_a", "3_aminoisobutanoate_b", "alanine_b")))]
### Remove missing compounds that missing in Plate 1 ODD batch...


p3e <- read.delim("mrmkit_output/2025-08-06-18-19-59 Jiangwen0513_2/quanttables/Plate3_EVEN_quant_table.txt", 
                  as.is = T, check.names = F)
p3e$batch <- 5
### Remove missing compounds that missing in Plate 1 ODD batch...
p3e <- p3e[-c(which(colnames(p3e) %in% c("3_aminoisobutanoate_a", "3_aminoisobutanoate_b", "alanine_b")))]
### Remove missing compounds that missing in Plate 1 ODD batch...


p3o <- read.delim("mrmkit_output/2025-08-06-18-19-59 Jiangwen0513_2/quanttables/Plate3_ODD_quant_table.txt",
                  as.is = T, check.names = F)
p3o$batch <- 6
### Remove missing compounds that missing in Plate 1 ODD batch...
p3o <- p3o[-c(which(colnames(p3o) %in% c("3_aminoisobutanoate_a", "3_aminoisobutanoate_b", "alanine_b")))]
### Remove missing compounds that missing in Plate 1 ODD batch...


### Read in S/N from raw data
setwd("mrmkit_output/2025-08-06-16-23-15 Jiangwen0513_5/quanttables/")

# snfile <- c("Plate1_EVEN_quant_table.txt", "Plate1_ODD_quant_table.txt",
#            "Plate2_EVEN_quant_table.txt", "Plate2_ODD_quant_table.txt",
#            "Plate3_EVEN_quant_table.txt", "Plate3_ODD_quant_table.txt")
# ncol <- ncol(read.delim("Plate1_ODD_quant_table.txt", as.is = T, check.names = F))
# 
# sn <- data.frame(matrix(ncol = ncol, nrow = 0))
# for (file in snfile){
#   sn <- rbind(sn, read.delim(file, as.is = T, check.names = F)[10,])
# }
# rm(ncol, snfile, file)

### Now that ODD quant tables have missing compounds, every quant table need to be read separately...
sn1e <- read.delim("Plate1_EVEN_quant_table.txt", as.is = T, check.names = F)[10, ]
sn1e$batch <- 1
sn1e <- sn1e[-c(which(colnames(sn1e) %in% c("3_aminoisobutanoate_a", "3_aminoisobutanoate_b", "alanine_b")))]

sn1o <- read.delim("Plate1_ODD_quant_table.txt", as.is = T, check.names = F)[10, ]
sn1o$batch <- 2

sn2e <- read.delim("Plate2_EVEN_quant_table.txt", as.is = T, check.names = F)[10, ]
sn2e$batch <- 3
sn2e <- sn2e[-c(which(colnames(sn2e) %in% c("3_aminoisobutanoate_a", "3_aminoisobutanoate_b", "alanine_b")))]

sn2o <- read.delim("Plate2_ODD_quant_table.txt", as.is = T, check.names = F)[10, ]
sn2o$batch <- 4
sn2o <- sn2o[-c(which(colnames(sn2o) %in% c("3_aminoisobutanoate_a", "3_aminoisobutanoate_b", "alanine_b")))]

sn3e <- read.delim("Plate3_EVEN_quant_table.txt", as.is = T, check.names = F)[10, ]
sn3e$batch <- 5
sn3e <- sn3e[-c(which(colnames(sn3e) %in% c("3_aminoisobutanoate_a", "3_aminoisobutanoate_b", "alanine_b")))]

sn3o <- read.delim("Plate3_ODD_quant_table.txt", as.is = T, check.names = F)[10, ]
sn3o$batch <- 6
sn3o <- sn3o[-c(which(colnames(sn3o) %in% c("3_aminoisobutanoate_a", "3_aminoisobutanoate_b", "alanine_b")))]
### Now that ODD quant tables have missing compounds, every quant table need to be read separately...

sn <- rbind(sn1e, sn1o, sn2e, sn2o, sn3e, sn3o)
# sn$batch <- c(1, 2, 3, 4, 5, 6)
sn <- sn[-3]

setwd(".")  # path made relative for public copy


### calculate median S/N with the 6 S/N values
msn <- c("Signal-to-noise(6median)","NA")
for(i in 3:ncol(sn)){
  msn <- c(msn, median(sn[[i]]))
}
rm(i)

sn <- rbind(sn, msn)
#** remember to remove those ISTDs cos (0,0) include ISTDS ####
rm(msn)

# write.csv(sn, file = "otpt/SignalToNoise_ratio_11.csv", row.names = F)

# Now filter out using S/N, keep those with S/N >= 5
sn_pass <- colnames(sn[7, ][-(1:2)])[which(as.numeric(sn[7, ][-(1:2)]) >= 5)]
### ^^^ inter-batch correction won't affect S/N of specific batches, so these S/N won't change after inter-batch corretcion.


### manual inter-batch correction - grand median
quant_bf <- rbind(p1e[11:nrow(p1e),], p1o[11:nrow(p1o),],
                  p2e[11:nrow(p2e),], p2o[11:nrow(p2o),],
                  p3e[11:nrow(p3e),], p3o[11:nrow(p3o),])
#^^this is the one whose intra-batch corection was done 
#by MRMKit automaticaly


# quant_bf <- read.delim("mrmkit_output/to_be_rerun", 
#                        as.is = T, check.names = F)
# quant_bf <- cbind(quant_bf[1:2], quant_bf[4], quant_bf[9:ncol(quant_bf)])
# colnames(quant_bf)[[3]] <- "type"
#^^this is the one whose intra-batch correction was done 
#by a labmate, manually using kernel regression


#** calculate grand median
gm <- NULL
for(i in 4:ncol(quant_bf)){
  gm <- c(gm, median(quant_bf[[i]], na.rm = T))
}
rm(i)
gm <- cbind(colnames(quant_bf)[4:ncol(quant_bf)], gm)

#** batch-correction using grand median
quant_af <- NULL
for (b in unique(quant_bf$batch)){ #* b: per batch
  batch_set <- subset(quant_bf, batch == b)
  for (a in 4:ncol(batch_set)){ #* a: per analyte
    batch_median <- median(batch_set[[a]], na.rm = T)
    batch_set[[a]] <- batch_set[[a]] * as.numeric(gm[[a-3, 2]]) / batch_median
  }
  rm(a, batch_median)
  
  quant_af <- rbind(quant_af, batch_set)
  #** quant_af is the quant table after batch correction
}
rm(b, batch_set)
### visualization: before and after batch correction comparison

# write.csv(quant_bf, file = "otpt/quant_table_before_InterBatch_MRMKit.csv", row.names = F)
# write.csv(quant_af, file = "otpt/quant_table_after_InterBatch_MRMKit.csv", row.names = F)


### cov, d-ratio and S/N filtering
af_berry <- quant_af %>% group_by(batch) %>% group_split

#* create qc-metric to contain cov and d-ratio values
CoV_metric <- D_ratio_metric <- data.frame(matrix(ncol = ncol(quant_af)-3, nrow = 0))
colnames(CoV_metric) <- colnames(D_ratio_metric) <- colnames(quant_af)[-(1:3)]

### calculate cov = 100 * SD(BQC)/Mean(BQC) AND d-ratio = 100 * SD(BQC)/SD(Samples)
#*keep manual calculation consistent with MRMKit: calculate for each batch separately, then take the median value
for(k in 1:length(af_berry)){

crt_batch <- unique(af_berry[[k]]$batch)
sub_bqc <- af_berry[[k]][which(af_berry[[k]]$type == "BQC"), ]
sub_sample <- af_berry[[k]][which(af_berry[[k]]$type == "SAMPLE"), ]

CoV <- data.frame(matrix(ncol = 0, nrow = 1))
D_ratio <- data.frame(matrix(ncol = 0, nrow = 1))

for(n in 4:ncol(sub_bqc)){
  crt_name <- colnames(sub_bqc)[[n]]
  
  bqc_value <- as.numeric(unlist(sub_bqc[, n]))
  sample_value <- as.numeric(unlist(sub_sample[, n]))
  
  crt_cov <- 100 * sd(bqc_value) / mean(bqc_value)
  crt_dratio <- 100 * sd(bqc_value) / sd(sample_value)
  
  
  CoV <- cbind(CoV, crt_cov)
  D_ratio <- cbind(D_ratio, crt_dratio)
  colnames(CoV)[[ncol(CoV)]] <- colnames(D_ratio)[[ncol(D_ratio)]] <- crt_name
}
rm(n, crt_name, bqc_value, sample_value, crt_cov, crt_dratio)

CoV_metric <- rbind(CoV_metric, CoV)
rownames(CoV_metric)[[(nrow(CoV_metric))]] <- paste0("CoV", " (Batch ", crt_batch, ")")
D_ratio_metric <- rbind(D_ratio_metric, D_ratio)
rownames(D_ratio_metric)[[nrow(D_ratio_metric)]] <- paste0("D_ratio", " (Batch ", crt_batch, ")")

}

rm(k, crt_batch, sub_bqc, sub_sample)

###Calculate median CoV and D-ratio
CoV_metric <- rbind(CoV_metric, apply(CoV_metric, 2, median, na.rm = TRUE))
rownames(CoV_metric)[[nrow(CoV_metric)]] <- "CoV (Median)"

D_ratio_metric <- rbind(D_ratio_metric, apply(D_ratio_metric, 2, median, na.rm = TRUE))
rownames(D_ratio_metric)[[nrow(D_ratio_metric)]] <- "D_ratio (Median)"

#** QC metrics (cov and d-ratio) filtering, keep those with CoV <= 30 OR D-ratio <= 50
cov_pass <- colnames(CoV_metric)[which(CoV_metric[nrow(CoV_metric), ] <= 30)]
dratio_pass <- colnames(D_ratio_metric)[which(D_ratio_metric[nrow(D_ratio_metric), ] <= 50)]

analytes_pass <- intersect(union(cov_pass, dratio_pass), sn_pass)

quant_af <- quant_af[, c(1:3, which(colnames(quant_af) %in% analytes_pass))]
