library(tidyverse)
library(data.table)

#dna barcode level data
dna_A <- fread("./data/barcode_level_data/SE_subsampled/EoE_repA_analysis_w_tag.txt")
dna_B <- fread("./data/barcode_level_data/SE_subsampled/EoE_repB_analysis_w_tag.txt")
dna_C <- fread("./data/barcode_level_data/SE_subsampled/EoE_repC_analysis_w_tag.txt")
dna_D <- fread("./data/barcode_level_data/SE_subsampled/EoE_repD_analysis_w_tag.txt")
dna_E <- fread("./data/barcode_level_data/SE_subsampled/EoE_repE_analysis_w_tag.txt")


#rna barcode level data TE7
rna_A <- fread("./data/barcode_level_data/MPRA-1_analysis_w_tag.txt")
rna_B <- fread("./data/barcode_level_data/MPRA-2_analysis_w_tag.txt")
rna_C <- fread("./data/barcode_level_data/MPRA-3_analysis_w_tag.txt")
rna_D <- fread("./data/barcode_level_data/MPRA-4_analysis_w_tag.txt")
rna_E <- fread("./data/barcode_level_data/MPRA-5_analysis_w_tag.txt")



#combine all the barcode and unique for a list
barcode <- bind_rows(
dna_A[,1:2], dna_B[,1:2], dna_C[,1:2], dna_D[,1:2],dna_E[,1:2],
rna_A[,1:2], rna_B[,1:2], rna_C[,1:2], rna_D[,1:2], rna_E[,1:2]
) %>% unique()


barcode_rank <- barcode %>%
arrange(reference_id) %>%
group_by(reference_id) %>%
mutate(barcode=row_number())


#dna_combine
dna <- dna_A[,2:3] %>%
full_join(dna_B[,2:3], by = "tag sequence", suffix= c(".A", ".B")) %>%
full_join(dna_C[,2:3], by = "tag sequence", suffix= c("", ".C")) %>%
full_join(dna_D[,2:3], by = "tag sequence", suffix= c("", ".D")) %>%
full_join(dna_E[,2:3], by = "tag sequence", suffix= c("", ".E")) %>%
rename(
	sampleA = "SE count.A", 
	sampleB="SE count.B", 
	sampleC="SE count", 
	sampleD="SE count.D", 
	sampleE="SE count.E",
)

#rna_combine
rna <- rna_A[,2:3] %>%
full_join(rna_B[,2:3], by = "tag sequence", suffix= c(".A", ".B")) %>%
full_join(rna_C[,2:3], by = "tag sequence",suffix= c("", ".C")) %>%
full_join(rna_D[,2:3], by = "tag sequence", suffix= c("", ".D")) %>%
full_join(rna_E[,2:3], by = "tag sequence", suffix= c("", ".E")) %>%
rename(
	sampleA = "SE count.A", 
	sampleB= "SE count.B", 
	sampleC= "SE count", 
	sampleD= "SE count.D", 
	sampleE= "SE count.E",
)


#add barcode and remove the tag sequence
dna_barcode <- full_join(dna, barcode_rank, by = "tag sequence") %>% select(-`tag sequence`)
rna_barcode <- full_join(rna, barcode_rank, by = "tag sequence") %>% select(-`tag sequence`)

#writeRDS
write_rds(dna_barcode, "./data/R_dataset/dna_barcode_level.RDS")
write_rds(rna_barcode, "./data/R_dataset/rna_barcode_level.RDS")
