library(fasthplus)
library(here)
library(usethis)
library(purrr)
library(ggplot2)
library(dplyr)
library(reshape2)

#call external scripts
source(file.path(here::here("simulations","01_simulate_data.R")))

## output directories
dir_plots <- here::here("simulation_figures", "02_hpe_hpb_diff" )
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)