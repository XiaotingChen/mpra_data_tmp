#setwd("/data/weirauchlab/team/ches2d/OngoingProjects/EoE_MPRA_paper_Molly/__revision_072023__/TE7_vs_control")

library(MPRAnalyze)
library(tidyverse)
library(BiocParallel)


dna_data <- readRDS("./data/R_dataset/dna_data_allelic.RDS")
rna_data <- readRDS("./data/R_dataset/rna_data_allelic.RDS")
colanno_data <- readRDS("./data/R_dataset/colanno_data_allelic.RDS")

colanno_data$barcode <- as.factor(colanno_data$barcode)
colanno_data$condition <- as.factor(colanno_data$type)
colanno_data$batch <- as.factor(colanno_data$batch)

# par.backend <- bpparam()
par.backend <- MulticoreParam(workers=16)

obj <- MpraObject(
	dnaCounts = dna_data,
	rnaCounts = rna_data,
    dnaAnnot = colanno_data,
	rnaAnnot = colanno_data,
    controls = NA_integer_ ,
	BPPARAM = par.backend
)


#compare two groups
obj <- analyzeComparative(obj = obj,
                              dnaDesign = ~ condition ,
                              rnaDesign = ~ condition,
                              reducedDesign = ~ 1
                              )



res <- testLrt(obj)

# res.epi<-res.epi %>%
# mutate(as.data.frame(rownames(res.epi))) %>%
# rename(label = "rownames(res.epi)") %>%
# separate(label, c("rsid", "type","allele"), remove = T, sep ="_")

saveRDS(obj, "./results/allelic_comparative_obj.RDS")
write.table(res,"./results/allelic_comparative.csv",sep="\t")
#write.csv(obj.res, "./results/EoE_TE7_vs_control_comparison_scale.csv")
