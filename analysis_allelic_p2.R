library(MPRAnalyze)
library(tidyverse)
library(BiocParallel)


# load dummy coded data
dna_eoe=read.csv(file="./data/R_dataset/dummy_coded_dna_eoe.csv",sep=" ")
rna_eoe=read.csv(file="./data/R_dataset/dummy_coded_rna_eoe.csv",sep=" ")

# sum count filter plus new columns
dna_eoe_clean <- dna_eoe %>%
replace(., is.na(.), 0) %>%
mutate(
	count_sum = sampleA + sampleB + sampleC + sampleD + sampleE,
	reference_id_barcode = paste(reference_id, barcode ,sep = "+"),
) %>%
filter(count_sum >=1) %>%
group_by(reference_id) %>%
mutate(new_barcode=row_number()) %>%
separate(reference_id, c("rsid", "type","allele"), remove = F, sep ="_") %>%
mutate(new_barcode_extended = paste(new_barcode, type, allele, sep="+"))


# > dna_eoe_clean[1:5,]
# # A tibble: 5 × 14
# # Groups:   reference_id [1]
#   sampleA sampleB sampleC sampleD sampleE reference_id        rsid  type  allele
#     <dbl>   <dbl>   <dbl>   <dbl>   <dbl> <chr>               <chr> <chr> <chr>
# 1       1       1       1       1       1 rs10038058_Non-Ref… rs10… Non-… G
# 2       2       2       2       2       2 rs10038058_Non-Ref… rs10… Non-… G
# 3       3       3       3       3       3 rs10038058_Non-Ref… rs10… Non-… G
# 4       4       4       4       4       4 rs10038058_Non-Ref… rs10… Non-… G
# 5       1       1       1       1       1 rs10038058_Non-Ref… rs10… Non-… G
# # ℹ 5 more variables: barcode <int>, count_sum <dbl>,
# #   reference_id_barcode <chr>, new_barcode <int>, new_barcode_w_type <chr>

# > colnames(dna_eoe_clean)
#  [1] "sampleA"              "sampleB"              "sampleC"
#  [4] "sampleD"              "sampleE"              "reference_id"
#  [7] "rsid"                 "type"                 "allele"
# [10] "barcode"              "count_sum"            "reference_id_barcode"
# [13] "new_barcode"          "new_barcode_extended"

# old arrange function
#arrange the data for MPRAanalyze
# arrange_fun <- function(x, y) {
#   colnames(x) <- c("sample", "reference_id", "barcode" )
#   wider_format <- x %>% pivot_wider(names_from = barcode, values_from = sample)
#   colnames(wider_format) <- paste(colnames(wider_format), y, sep = "_")
#   return(wider_format)
# } #give sample, reference_id and barcode
#
arrange_fun <- function(x, y) {
  colnames(x) <- c("sample", "rsid", "new_barcode_extended" )
  wider_format <- x %>% pivot_wider(names_from = new_barcode_extended, values_from = sample)
  colnames(wider_format) <- paste(colnames(wider_format), y, sep = "+")
  return(wider_format)
} #give sample, reference_id and barcode


#dna_arrange | sample:1-5, rsid:7, new_barcode_extended:14
dna_sampleA <- arrange_fun(dna_eoe_clean[,c(1,7,14)], "NA_A") %>% replace(., is.na(.), 0)
dna_sampleB <- arrange_fun(dna_eoe_clean[,c(2,7,14)], "NA_B") %>% replace(., is.na(.), 0)
dna_sampleC <- arrange_fun(dna_eoe_clean[,c(3,7,14)], "NA_C") %>% replace(., is.na(.), 0)
dna_sampleD <- arrange_fun(dna_eoe_clean[,c(4,7,14)], "NA_D") %>% replace(., is.na(.), 0)
dna_sampleE <- arrange_fun(dna_eoe_clean[,c(5,7,14)], "NA_E") %>% replace(., is.na(.), 0)

# > dim(dna_sampleA)
# [1]   551 18832
# > dna_sampleA[1:5,1:5]
# # A tibble: 5 × 5
#   `rsid+NA_A` `1+Non-Ref+G+NA_A` `2+Non-Ref+G+NA_A` `3+Non-Ref+G+NA_A`
#   <chr>                    <dbl>              <dbl>              <dbl>
# 1 rs10038058                   1                  2                  3
# 2 rs10038177                   0                  0                  0
# 3 rs10043631                   0                  0                  0
# 4 rs10045255                   1                  2                  1
# 5 rs10050834                   0                  0                  0



dna_arrange <- left_join(dna_sampleA, dna_sampleB, c("rsid+NA_A" = "rsid+NA_B")) %>%
left_join(dna_sampleC, c("rsid+NA_A" = "rsid+NA_C")) %>%
left_join(dna_sampleD, c("rsid+NA_A" = "rsid+NA_D")) %>%
left_join(dna_sampleE, c("rsid+NA_A" = "rsid+NA_E"))

dna_data <- dna_arrange[,-1] %>% as.matrix()
rownames(dna_data) <- dna_arrange$"rsid+NA_A"


# rna

rna_eoe_clean <- rna_eoe %>%
replace(., is.na(.), 0) %>%
mutate(
	reference_id_barcode = paste(reference_id, barcode ,sep = "+")
) %>%
filter(reference_id_barcode %in% dna_eoe_clean$reference_id_barcode) %>%
left_join(dna_eoe_clean[, c("reference_id_barcode", "rsid","new_barcode_extended")],
"reference_id_barcode")

# > dim(rna_eoe_clean)
# [1] 556040     10
# > colnames(rna_eoe_clean)
#  [1] "sampleA"              "sampleB"              "sampleC"
#  [4] "sampleD"              "sampleE"              "reference_id"
#  [7] "barcode"              "reference_id_barcode" "rsid"
# [10] "new_barcode_extended"


#rna_arrange
rna_sampleA <- arrange_fun(rna_eoe_clean[,c(1,9,10)], "NA_A") %>% replace(., is.na(.), 0)
rna_sampleB <- arrange_fun(rna_eoe_clean[,c(2,9,10)], "NA_B") %>% replace(., is.na(.), 0)
rna_sampleC <- arrange_fun(rna_eoe_clean[,c(3,9,10)], "NA_C") %>% replace(., is.na(.), 0)
rna_sampleD <- arrange_fun(rna_eoe_clean[,c(4,9,10)], "NA_D") %>% replace(., is.na(.), 0)
rna_sampleE <- arrange_fun(rna_eoe_clean[,c(5,9,10)], "NA_E") %>% replace(., is.na(.), 0)


rna_arrange <- left_join(rna_sampleA, rna_sampleB, c("rsid+NA_A" = "rsid+NA_B")) %>%
left_join(rna_sampleC, c("rsid+NA_A" = "rsid+NA_C")) %>%
left_join(rna_sampleD, c("rsid+NA_A" = "rsid+NA_D")) %>%
left_join(rna_sampleE, c("rsid+NA_A" = "rsid+NA_E"))

rna_data <- rna_arrange[,-1] %>% as.matrix()
rownames(rna_data) <- rna_arrange$"rsid+NA_A"


# check if all column names are the same
# all(colnames(dna_arrange) == colnames(rna_arrange))

colanno <- as.data.frame(colnames(rna_data)) %>%
rename(barcode_info = "colnames(rna_data)") %>%
separate(barcode_info, c("barcode", "type","allele","batch"), remove = F, sep="\\+")



colanno_data <- colanno[,-1]
colanno_data$barcode <- as.factor(colanno_data$barcode)
colanno_data$type <- as.factor(colanno_data$type)
colanno_data$batch <- as.factor(colanno_data$batch)
rownames(colanno_data) <- colanno$barcode_info

# output
write_rds(dna_data, "./data/R_dataset/dna_data_allelic.RDS")
write_rds(rna_data, "./data/R_dataset/rna_data_allelic.RDS")
write_rds(colanno_data, "./data/R_dataset/colanno_data_allelic.RDS")

