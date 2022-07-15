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
wp_raw <- read_excel('FRZB_Dataset.xlsx', sheet = 'Phospho', na = c("", "NA")) %>%
  select(c(`Accession`,
           `Abundances (Grouped)`)) %>%
  na.omit()