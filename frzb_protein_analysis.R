# About
# This script reads in raw data > preprocesses data > to form a data matrix
# Data matrix is used for further stats analysis

# ============== Clears Memory ======
# clears all objects includes hidden objects
rm(list = ls(all.names = TRUE)) 

# frees up memory and reports the memory usage.
gc() 

# ============== Loads Packages =======
library(readxl)
library(dplyr)
library(data.table)
library(berryFunctions)
library(destiny)
library(tidyr)
library(stringr)
library(ggplot2)
library(janitor)
library(IMIFA)
library(tidyverse)

# ============== 1. Read & Selects from Excel File ======

# reads in raw data
wp_raw <- read_excel('FRZB_Dataset.xlsx', sheet = 'Whole', na = c("", "NA")) %>%
  select(c(`Accession`,
           `Abundances (Grouped)`)) %>%
  na.omit()

# ============== 2. Manipulates Data  ===============
# manipulates data
{
  # splits string into columns
  wp_raw_split <- str_split_fixed(as.character(wp_raw$`Abundances (Grouped)`), ';',16)
  
  # adds columns to original data
  wp_raw <- cbind(wp_raw, wp_raw_split)
  colnames(wp_raw) <- c("Accession",
                        "Abundances (Grouped)",
                        "S1F",
                        'S1V',
                        'S2F',
                        'S2V',
                        'S3F',
                        'S3V')
  
  # removes abundance column
  wp_clean <- wp_raw[, -2]
  wp_clean <- wp_clean %>% discard(~all(is.na(.) | . ==""))
}

# exports abundance matrix to csv
fwrite(wp_clean, "Whole_Protein_Accession.csv", sep = ",")

# ============== 3. Combines Uniprot Data To Combined Matrix =====
# reads in Gene Symbol table downloaded from Uniprot
gene_symbol <- fread("Whole_Protein_Accession_Map.tsv")

# renames column headers
colnames(gene_symbol) <- c("From", "Accession", "Gene Symbol") 

# removes first column
gene_symbol <- gene_symbol[,-1]

# splits gene symbol to return only the first 
gene_symbol$`Gene Symbol` <- sapply(strsplit(gene_symbol$`Gene Symbol`," "), `[`, 1)

# merges gene symbol column to main df
wp_clean_gs <- wp_clean %>%
  
  # merges gene symbol column to main df
  left_join(gene_symbol,
            by = "Accession") %>%
  
  distinct(`Accession`, .keep_all = TRUE) %>%

  
  # relocates columns and removes NAs
  relocate(c(`Accession`,
             `Gene Symbol`),
           .before = `S1F`) %>%
  na.omit() %>%
  
  # adds number to the end of duplicate gene symbols (ie Sptbn1-2)
  group_by(`Gene Symbol`) %>%
  mutate(`GS_count` = 1:n()) %>%
  mutate(`Gene Symbol` = ifelse(`GS_count` == 1, 
                                `Gene Symbol`, 
                                paste0(`Gene Symbol`, "-", `GS_count`))) %>%
  
  # removes unused columns
  select(-`GS_count`) 

# ============== 4. Exports Abundance Matrix =====

# exports grouped abundance matrix to csv
fwrite(wp_clean_gs, "wp_grouped.csv", sep = ",")

