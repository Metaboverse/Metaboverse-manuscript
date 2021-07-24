# Data sources:
# "rna_cancer_sample.tsv": Human Protein Atlas (https://www.proteinatlas.org/download/rna_cancer_sample.tsv.zip)
# "clinical.tsv": TCGA (https://portal.gdc.cancer.gov/projects/TCGA-LUAD; click Clinical download button)
library(survival)
library(survminer)


setwd("C:/Users/jorda/Desktop/projects/manuscript/rna_cancer_data/")

rpkm_data <- read.csv(file = 'rna_cancer_sample.tsv', sep='\t')
# Extract LUAD data
rpkm_data <- rpkm_data[rpkm_data$Cancer == 'LUAD',]
# Clean up Sample IDs
rpkm_data$Sample <- sapply(strsplit(rpkm_data$Sample, "-01A"), "[", 1)
gene_list <- unique(rpkm_data$Gene)
rpkm_data$FPKM <- as.numeric(rpkm_data$FPKM)

clinical_data <- read.csv(file = "clinical.tsv", sep='\t')

# Censor alive/dead
censor_data <- function(x) if (x == "Alive") {0} else {1}
clinical_data$status <- mapply(censor_data, clinical_data$vital_status)
clinical_data$status <- as.numeric(clinical_data$status)
clinical_data$days_to_death <- as.numeric(clinical_data$days_to_death)

# Censor death date for alive patients (N/A) with no final status
#my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
#this_max <- my.max(clinical_data$days_to_death)
#clinical_data$days_to_death[is.na(clinical_data$days_to_death)] <- this_max

# days_to_last_follow_up
clinical_data$days_to_death <- ifelse(
  is.na(clinical_data$days_to_death), 
  clinical_data$days_to_last_follow_up, 
  clinical_data$days_to_death)

clinical_sample <- clinical_data[ , c("case_submitter_id", "ajcc_pathologic_stage", "status", "days_to_death", "days_to_last_follow_up")] 
clinical_sample <- clinical_sample[clinical_sample$days_to_death > 10 & clinical_sample$days_to_last_follow_up != "'--", ]
clinical_sample <- clinical_sample[clinical_sample$days_to_death != "'--", ]

clinical_sample$status <- as.numeric(as.character(clinical_sample$status))
clinical_sample$days_to_death <- as.numeric(as.character(clinical_sample$days_to_death))
clinical_sample$days_to_last_follow_up <- as.numeric(as.character(clinical_sample$days_to_last_follow_up))

clinical_sample$stage <- clinical_sample$ajcc_pathologic_stage
clinical_sample <- clinical_sample[clinical_sample$stage != "'--", ]
clinical_sample$stage[
  clinical_sample$stage == "Stage I" | 
  clinical_sample$stage == "Stage IA" | 
  clinical_sample$stage == "Stage IB"] <- 1
clinical_sample$stage[
  clinical_sample$stage == "Stage II" | 
  clinical_sample$stage == "Stage IIA" | 
  clinical_sample$stage == "Stage IIB"] <- 2
clinical_sample$stage[
  clinical_sample$stage == "Stage III" | 
  clinical_sample$stage == "Stage IIIA" | 
  clinical_sample$stage == "Stage IIIB"] <- 3
clinical_sample$stage[
  clinical_sample$stage == "Stage IV" | 
  clinical_sample$stage == "Stage IVA" | 
  clinical_sample$stage == "Stage IVB"] <- 4
  

clinical_sample_early <- clinical_sample[clinical_sample$ajcc_pathologic_stage == "Stage I" | clinical_sample$ajcc_pathologic_stage == "Stage IA" | clinical_sample$ajcc_pathologic_stage == "Stage IB", ]
clinical_sample_late <- clinical_sample[clinical_sample$ajcc_pathologic_stage != "Stage I" & clinical_sample$ajcc_pathologic_stage != "Stage IA" & clinical_sample$ajcc_pathologic_stage != "Stage IB", ]


# Extract LUAD data for SMS (ENSG00000102172) (loop this)
median_pvals <-list()
optimized_pvals <- list()

for (this_gene in gene_list) {
  tryCatch({
    data_gene <- rpkm_data[rpkm_data$Gene == this_gene,]
    # Merged mixed data tables
    merged_data_gene <- merge(
      data_gene,
      clinical_sample,
      by.x = "Sample",
      by.y = "case_submitter_id")

    this_median <- median(merged_data_gene$FPKM)
    merged_data_gene$gene_status <- ifelse(merged_data_gene$FPKM > this_median, "high", "low")

    km_trt_fit <- survfit(Surv(days_to_death, status) ~ gene_status, data=merged_data_gene)

    opt_cutoff <- surv_cutpoint(merged_data_gene, time = "days_to_death", event = "status", variables = c("FPKM"))
    # print(c(this_gene, opt_cutoff))
    opt_categories <- surv_categorize(opt_cutoff)
    opt_fit <- survfit(Surv(days_to_death, status) ~ FPKM, data = opt_categories)

    median_pvals[this_gene] <- surv_pvalue(km_trt_fit)$pval
    optimized_pvals[this_gene] <- surv_pvalue(opt_fit)$pval
  }, error = function (condition) {
    print(this_gene)
    #print(condition)
    #break
  })
}






# this_gene <- "ENSG00000102172" # SMS
this_gene <- "ENSG00000102172"
data_gene <- rpkm_data[rpkm_data$Gene == this_gene,]
# Merged mixed data tables
merged_data_gene <- merge(
  data_gene,
  clinical_sample,
  by.x = "Sample",
  by.y = "case_submitter_id")

this_median <- median(merged_data_gene$FPKM)
merged_data_gene$gene_status <- ifelse(merged_data_gene$FPKM > this_median, "high", "low")

km_trt_fit <- survfit(Surv(days_to_death, status) ~ gene_status, data=merged_data_gene)

opt_cutoff <- surv_cutpoint(merged_data_gene, time = "days_to_death", event = "status", variables = c("FPKM"))
opt_categories <- surv_categorize(opt_cutoff)
opt_fit <- survfit(Surv(days_to_death, status) ~ FPKM, data = opt_categories)

median_pvals[this_gene] <- surv_pvalue(km_trt_fit)$pval
optimized_pvals[this_gene] <- surv_pvalue(opt_fit)$pval
sms_val <- surv_pvalue(opt_fit)$pval

ggsurvplot(
  opt_fit,
  risk.table = TRUE,
  conf.int=TRUE,
  pval=sms_val,
  pval.method=TRUE,
  palette=c("dodgerblue2", "orchid2"),
  title="Kaplan-Meier Curve for Lung Cancer Survival",
  risk.table.height=.25)

print(this_gene)
res.cox <- coxph(Surv(days_to_death, status) ~ gene_status + ajcc_pathologic_stage, data=merged_data_gene)
summary(res.cox)




# this_gene <- "ENSG00000168237" # GLYCTK
this_gene <- "ENSG00000168237"
data_gene <- rpkm_data[rpkm_data$Gene == this_gene,]
# Merged mixed data tables
merged_data_gene <- merge(
  data_gene,
  clinical_sample,
  by.x = "Sample",
  by.y = "case_submitter_id")

this_median <- median(merged_data_gene$FPKM)
merged_data_gene$gene_status <- ifelse(merged_data_gene$FPKM > this_median, "high", "low")

km_trt_fit <- survfit(Surv(days_to_death, status) ~ gene_status, data=merged_data_gene)

opt_cutoff <- surv_cutpoint(merged_data_gene, time = "days_to_death", event = "status", variables = c("FPKM"))
# print(c(this_gene, opt_cutoff))
opt_categories <- surv_categorize(opt_cutoff)
opt_fit <- survfit(Surv(days_to_death, status) ~ FPKM, data = opt_categories)

median_pvals[this_gene] <- surv_pvalue(km_trt_fit)$pval
optimized_pvals[this_gene] <- surv_pvalue(opt_fit)$pval
glyctk_val <- surv_pvalue(opt_fit)$pval

ggsurvplot(
  opt_fit,
  risk.table = TRUE,
  conf.int=TRUE,
  pval=glyctk_val,
  pval.method=TRUE,
  palette=c("dodgerblue2", "orchid2"),
  title="Kaplan-Meier Curve for Lung Cancer Survival",
  risk.table.height=.25)

print(this_gene)
res.cox <- coxph(Surv(days_to_death, status) ~ gene_status, data=merged_data_gene)
summary(res.cox)




# Plot distribution
bump <- 1e-5
hist(
  -1 * log10(as.numeric(optimized_pvals) + bump),
  breaks=150,
  main="LUAD optimized survival distribution",
  xlab=expression("-log"[10]*"(log-rank p-value)"))
abline(v = -1 * log10(sms_val + bump), col="#1b9e77", lwd=3, lty=2)
abline(v = -1 * log10(glyctk_val + bump), col="#7570b3", lwd=3, lty=2)



# Get ranked percentile
perc.rank <- function(x, xo)  length(x[x <= xo])/length(x)*100
optimized_pvals_sorted <- optimized_pvals[order(unlist(optimized_pvals), decreasing=FALSE)]
perc.rank(optimized_pvals_sorted, sms_val) #0.01845855
perc.rank(optimized_pvals_sorted, glyctk_val) #0.1791883



optimized_pvals_bh <- p.adjust(
  optimized_pvals_sorted, 
  method = "BH", 
  n = length(optimized_pvals_sorted))

hist(
  -1 * log10(as.numeric(optimized_pvals_bh) + bump),
  breaks=150,
  main="LUAD optimized survival distribution (BH-corrected)",
  xlab=expression("-log"[10]*"(log-rank p-value)"))
abline(v = -1 * log10(optimized_pvals_bh['ENSG00000102172'] + bump), col="#1b9e77", lwd=3, lty=2)
abline(v = -1 * log10(optimized_pvals_bh['ENSG00000168237'] + bump), col="#7570b3", lwd=3, lty=2)



bump <- 1
optimized_pvals_bon <- p.adjust(
  optimized_pvals_sorted, 
  method = "bonferroni", 
  n = length(optimized_pvals_sorted))

hist(
  -1 * log10(as.numeric(optimized_pvals_bon) + bump),
  breaks=150,
  main="LUAD optimized survival distribution (Bonferroni-corrected)",
  xlab=expression("-log"[10]*"(log-rank p-value)"))
abline(v = -1 * log10(optimized_pvals_bon['ENSG00000102172'] + bump), col="#1b9e77", lwd=3, lty=2)
abline(v = -1 * log10(optimized_pvals_bon['ENSG00000102172'] + bump), col="#7570b3", lwd=3, lty=2)




sessionInfo()

session_output <- "
> sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19042)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] survminer_0.4.9 ggpubr_0.4.0    ggplot2_3.3.3   survival_3.2-11

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.6            mvtnorm_1.1-1         lattice_0.20-41       tidyr_1.1.2           zoo_1.8-8
 [6] digest_0.6.27         assertthat_0.2.1      utf8_1.2.1            R6_2.5.0              cellranger_1.1.0
[11] backports_1.2.1       pillar_1.6.0          rlang_0.4.10          curl_4.3              readxl_1.3.1
[16] data.table_1.13.6     car_3.0-10            Matrix_1.2-18         maxstat_0.7-25        labeling_0.4.2
[21] splines_4.0.3         stringr_1.4.0         foreign_0.8-80        munsell_0.5.0         gridtext_0.1.4
[26] tinytex_0.31          broom_0.7.6           compiler_4.0.3        xfun_0.22             pkgconfig_2.0.3
[31] ggtext_0.1.1          tidyselect_1.1.0      tibble_3.0.6          gridExtra_2.3         km.ci_0.5-2
[36] rio_0.5.26            fansi_0.4.2           crayon_1.4.1          dplyr_1.0.5           withr_2.4.2
[41] grid_4.0.3            xtable_1.8-4          gtable_0.3.0          lifecycle_1.0.0       DBI_1.1.1
[46] magrittr_2.0.1        KMsurv_0.1-5          scales_1.1.1          zip_2.1.1             stringi_1.5.3
[51] carData_3.0-4         farver_2.1.0          ggsignif_0.6.1        xml2_1.3.2            ellipsis_0.3.1
[56] survMisc_0.5.5        generics_0.1.0        vctrs_0.3.6           openxlsx_4.2.3        tools_4.0.3
[61] forcats_0.5.1         glue_1.4.2            markdown_1.1          purrr_0.3.4           hms_1.0.0
[66] abind_1.4-5           colorspace_2.0-0      rstatix_0.7.0         knitr_1.33            haven_2.3.1
[71] exactRankTests_0.8-32
Warning messages:
1: package 'survival' was built under R version 4.0.5
2: package 'survminer' was built under R version 4.0.5
3: In eval(ei, envir) : NAs introduced by coercion
"
