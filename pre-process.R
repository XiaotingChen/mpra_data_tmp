# package needed ----------------------------------------------------------

library(readxl)
#devtools::install_git("https://tfwebdev.research.cchmc.org/gitlab/lux2ht/mpraprofiler.git", ,credentials = git2r::cred_user_pass("lux2ht", getPass::getPass()))
library(mpraprofiler)
library("DESeq2")
library(tidyverse)
library(corrplot)
library(cowplot)

# raw sample counts -------------------------------------------------------
###### import re-sampling data
file_list <- list.files("./data/",pattern = '*.xlsx')
file_loc <- paste0("./data/", file_list)
df_list <- lapply(file_loc, read_xlsx)
mpra_name <- c("Ctrl","repA", "repB", "repC", "repD","repE")
names(df_list) <- mpra_name
df <- bind_rows(df_list, .id = "id")
df_use <- df %>% select(1,3,8) #
names(df_use) <- c("sample","oligo","counts")
count_table <- spread(df_use,sample, counts)

#check the correlation
# pairs(count_table_sle[,-1])
#cor(count_table[,-1])
#M <- floor(cor(count_table[,-1])*100)/100
#corrplot(M, method="color",  
#         type="lower",
#         addCoef.col = "black", # Add coefficient of correlation
#         tl.col="black", tl.srt = 0#Text label color and rotation
#)

#extract useful one
#read annotation table
seq_info <- read_xlsx("./data/annotation/Supplementary Data 1.xlsx", sheet = "MPRA library oligos")

#extract useful one
count_table_use <- count_table %>% filter(oligo %in% seq_info$Oligo_ID)

#save to others
saveRDS(count_table_use, "./data/R_dataset/Oligo_level_EoE_data.RDS")
