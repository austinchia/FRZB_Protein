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
pp_raw <- read_excel('FRZB_Dataset.xlsx', sheet = 'Phospho', na = c("", "NA")) %>%
  select(c(`Annotated Sequence`,
           `Modifications`,
           `Master Protein Accessions`,
           `Abundances (Grouped)`)) %>%
  na.omit()

# ============== 2. Manipulates Data  ===============
# manipulates data
{
  # splits string into columns
  pp_raw_split <- str_split_fixed(as.character(pp_raw$`Abundances (Grouped)`), ';',6)
  
  # adds columns to original data
  pp_raw <- cbind(pp_raw, pp_raw_split)
  colnames(pp_raw) <- c("Annotated Sequence",
                        "Modifications",
                        "Master Protein Accessions",
                        "Abundances (Grouped)",
                        "S1F",
                        'S1V',
                        'S2F',
                        'S2V',
                        'S3F',
                        'S3V')
  
  # removes abundance column
  pp_clean <- pp_raw[, -4]
}

# splits accession by ";" delimiter (ie "Q9JHU4-1; Q9JHU4" --> "Q9JHU4-1")
pp_clean$`Master Protein Accessions` <- sapply(strsplit(pp_clean$`Master Protein Accessions`,";"), `[`, 1)

# exports abundance matrix to csv
fwrite(pp_clean, "Phospho_Protein_Accession.csv", sep = ",")

# ============== 4. Combines Uniprot Data To Combined Matrix =====
# reads in Gene Symbol table downloaded from Uniprot
gene_symbol <- fread("Phospho_Protein_Accession_Map.tsv")

# renames column headers
colnames(gene_symbol) <- c("From", "Master Protein Accessions", "Gene Symbol") 

# splits gene symbol to return only the first 
gene_symbol$`Gene Symbol` <- sapply(strsplit(gene_symbol$`Gene Symbol`," "), `[`, 1)

# merges gene symbol column to main df
pp_clean_gs <- pp_clean %>%
  
  # merges gene symbol column to main df
  left_join(gene_symbol,
            by = "Master Protein Accessions") %>%
  
  # distinct(`Accession`, .keep_all = TRUE) %>%
  
  
  # relocates columns and removes NAs
  relocate(c(`Annotated Sequence`,
             `Modifications`,
             `Master Protein Accessions`,
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
  select(-c(`GS_count`,
         `From`)) 

# exports grouped abundance matrix to csv
fwrite(pp_clean_gs, "pp_grouped.csv", sep = ",")
