library(MPRAnalyze)
library(tidyverse)
library(BiocParallel)


#load data
#dna
dna <- readRDS("./data/R_dataset/dna_barcode_level.RDS") %>% arrange(reference_id,barcode)
rna <- readRDS("./data/R_dataset/rna_barcode_level.RDS") %>% arrange(reference_id, barcode)


#filter DNA to only EoE
count_table <- readRDS("./data/R_dataset/Oligo_level_EoE_data.RDS")
dna_eoe <- dna %>% filter(reference_id %in% count_table$oligo)
rna_eoe <- rna %>% filter(reference_id %in% count_table$oligo)


dna_eoe_clean <- dna_eoe %>%
replace(., is.na(.), 0) %>% 
mutate(
	count_sum = sampleA + sampleB + sampleC + sampleD + sampleE, 
	reference_id_barcode = paste(reference_id, barcode ,sep = "+")
) %>%
filter(count_sum >=1) %>%
group_by(reference_id) %>%
mutate(new_barcode=row_number())


rna_eoe_clean <- rna_eoe %>%
replace(., is.na(.), 0) %>%
mutate(
	reference_id_barcode = paste(reference_id, barcode ,sep = "+")
) %>%
filter(reference_id_barcode %in% dna_eoe_clean$reference_id_barcode) %>%
left_join(dna_eoe_clean[, c("reference_id_barcode", "new_barcode")], "reference_id_barcode")


# hist(dna_sle_clean$count_sum, breaks = 100, xlim = c(0,10))

#arrange the data for MPRAanalyze
arrange_fun <- function(x, y) {
  colnames(x) <- c("sample", "reference_id", "barcode" )
  wider_format <- x %>% pivot_wider(names_from = barcode, values_from = sample)
  colnames(wider_format) <- paste(colnames(wider_format), y, sep = "_")
  return(wider_format)
} #give sample, reference_id and barcode


#dna_arrange
dna_sampleA <- arrange_fun(dna_eoe_clean[,c(1,6,10)], "NA_A") %>% replace(., is.na(.), 0) 
dna_sampleB <- arrange_fun(dna_eoe_clean[,c(2,6,10)], "NA_B") %>% replace(., is.na(.), 0) 
dna_sampleC <- arrange_fun(dna_eoe_clean[,c(3,6,10)], "NA_C") %>% replace(., is.na(.), 0) 
dna_sampleD <- arrange_fun(dna_eoe_clean[,c(4,6,10)], "NA_D") %>% replace(., is.na(.), 0) 
dna_sampleE <- arrange_fun(dna_eoe_clean[,c(5,6,10)], "NA_E") %>% replace(., is.na(.), 0) 


dna_arrange <- left_join(dna_sampleA, dna_sampleB, c("reference_id_NA_A" = "reference_id_NA_B")) %>%
left_join(dna_sampleC, c("reference_id_NA_A" = "reference_id_NA_C")) %>%
left_join(dna_sampleD, c("reference_id_NA_A" = "reference_id_NA_D")) %>%
left_join(dna_sampleE, c("reference_id_NA_A" = "reference_id_NA_E"))

dna_data <- dna_arrange[,-1] %>% as.matrix()
rownames(dna_data) <- dna_arrange$reference_id_NA_A


#rna_arrange
rna_sampleA <- arrange_fun(rna_eoe_clean[,c(1,6,9)], "NA_A") %>% replace(., is.na(.), 0) 
rna_sampleB <- arrange_fun(rna_eoe_clean[,c(2,6,9)], "NA_B") %>% replace(., is.na(.), 0) 
rna_sampleC <- arrange_fun(rna_eoe_clean[,c(3,6,9)], "NA_C") %>% replace(., is.na(.), 0) 
rna_sampleD <- arrange_fun(rna_eoe_clean[,c(4,6,9)], "NA_D") %>% replace(., is.na(.), 0) 
rna_sampleE <- arrange_fun(rna_eoe_clean[,c(5,6,9)], "NA_E") %>% replace(., is.na(.), 0)


rna_arrange <- left_join(rna_sampleA, rna_sampleB, c("reference_id_NA_A" = "reference_id_NA_B")) %>%
left_join(rna_sampleC, c("reference_id_NA_A" = "reference_id_NA_C")) %>%
left_join(rna_sampleD, c("reference_id_NA_A" = "reference_id_NA_D")) %>%
left_join(rna_sampleE, c("reference_id_NA_A" = "reference_id_NA_E"))

rna_data <- rna_arrange[,-1] %>% as.matrix()
rownames(rna_data) <- rna_arrange$reference_id_NA_A


# check if all column names are the same
# all(colnames(dna_arrange) == colnames(rna_arrange))

colanno <- as.data.frame(colnames(rna_data)) %>%
rename(barcode_info = "colnames(rna_data)") %>%
separate(barcode_info, c("barcode", "condition","batch"), remove = F) %>%
mutate(condition) %>%
mutate(batch) 


colanno_data <- colanno[,-1]
colanno_data$barcode <- as.factor(colanno_data$barcode)
colanno_data$sample <- as.factor(colanno_data$condition)
colanno_data$batch <- as.factor(colanno_data$batch)
rownames(colanno_data) <- colanno$barcode_info

# output 
write_rds(dna_data, "./data/R_dataset/dna_data.RDS")
write_rds(rna_data, "./data/R_dataset/rna_data.RDS")
write_rds(colanno_data, "./data/R_dataset/colanno_data.RDS")

