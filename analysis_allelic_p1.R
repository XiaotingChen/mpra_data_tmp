library(MPRAnalyze)
library(tidyverse)
library(BiocParallel)


#load data
#dna
dna <- readRDS("./data/R_dataset/dna_barcode_level.RDS") %>% arrange(reference_id,barcode)
rna <- readRDS("./data/R_dataset/rna_barcode_level.RDS") %>% arrange(reference_id, barcode)

# > dim(dna)
# [1] 4850867       7
# > dna[1:5,]
#    sampleA sampleB sampleC sampleD sampleE        reference_id barcode
# 1:       3       3       3       3       3 rs1000309_Non-Ref_A       1
# 2:       2       2       2       2       2 rs1000309_Non-Ref_A       2
# 3:       3       3       3       3       3 rs1000309_Non-Ref_A       3
# 4:       3       3       3       3       3 rs1000309_Non-Ref_A       4
# 5:       5       5       5       5       5 rs1000309_Non-Ref_A       5


#filter DNA to only EoE oligos
count_table <- readRDS("./data/R_dataset/Oligo_level_EoE_data.RDS")
dna_eoe <- dna %>% filter(reference_id %in% count_table$oligo)
rna_eoe <- rna %>% filter(reference_id %in% count_table$oligo)

# export for adding dummy coding of labels
write.table(dna_eoe,file='./data/R_dataset/dna_eoe.csv')
write.table(rna_eoe,file='./data/R_dataset/rna_eoe.csv')

# call python code to get refined data
