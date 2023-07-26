setwd("/data/weirauchlab/team/ches2d/OngoingProjects/EoE_MPRA_paper_Molly/__revision_072023__/TE7_vs_control")

library(MPRAnalyze)
library(tidyverse)
library(BiocParallel)


dna_data <- readRDS("./data/R_dataset/dna_data.RDS")
rna_data <- readRDS("./data/R_dataset/rna_data.RDS")
colanno_data <- readRDS("./data/R_dataset/colanno_data.RDS")

colanno_data$barcode <- as.factor(colanno_data$barcode)
colanno_data$condition <- as.factor(colanno_data$condition)
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
obj <- analyzeQuantification(obj = obj,
                              dnaDesign = ~ 1 ,
                              rnaDesign = ~ 1)

alpha <- getAlpha(obj)
res.epi <- testEmpirical(obj = obj,
                         twoSided=TRUE
                         )

res.epi<-res.epi %>%
mutate(as.data.frame(rownames(res.epi))) %>%
rename(label = "rownames(res.epi)") %>%
separate(label, c("rsid", "type","allele"), remove = T, sep ="_")

saveRDS(obj, "./results/EoE_TE7_vs_control_obj.RDS")
write.table(res.epi,"./results/EoE_TE7_vs_control_quantification.csv",sep="\t")
#write.csv(obj.res, "./results/EoE_TE7_vs_control_comparison_scale.csv")
