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
wp_raw <- read_excel('NTU_Myopia_Dataset.xlsx', sheet = 'Whole_NTU', na = c("", "NA")) %>%
  select(c(`Accession`,
           `Abundances (Grouped)`)) %>%
  na.omit()

# ============== 2. Manipulates Data  ===============
# S1 Data manipulation
{
  # splits string into columns
  wp_raw_split <- str_split_fixed(as.character(wp_raw$`Abundances (Grouped)`), ';',16)
  
  # adds columns to original data
  wp_raw <- cbind(wp_raw, wp_raw_split)
  colnames(wp_raw) <- c("Accession", "Abundances (Grouped)",
                        "GC_S1",
                        'S1_LI_1hr','S1_LI_6hr','S1_LI_9hr','S1_LI_D1','S1_LI_D14','S1_LI_D3','S1_LI_D7','S1_NL_0hr','S1_NL_1hr','S1_NL_6hr','S1_NL_9hr','S1_NL_D1','S1_NL_D14','S1_NL_D3','S1_NL_D7')
  
  # removes abundance column
  wp_raw <- select(wp_raw, -c(`Abundances (Grouped)`)) %>%
    mutate(S1_LI_0hr = S1_NL_0hr) %>%
    relocate(S1_LI_0hr, .after = `Accession`) %>%
    relocate(S1_LI_D14, .after = `S1_LI_D7`) %>%
    relocate(S1_NL_D14, .after = `S1_NL_D7`)
  
}

