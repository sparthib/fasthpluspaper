#!/usr/bin/env Rscript

# install.packages('devtools')
# library(devtools)
# install_github(repo="ntdyjack/fasthplus", ref = "main")#load relevant packages
library(fasthplus)
library(here)
library(usethis)
library(purrr)
library(ggplot2)
library(dplyr)
library(reshape2)

##source scripts


## output directories
# dir_plots <- here::here("simulation_figures", "01_across_r" )
# dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

#code for simulating data 
#https://github.com/stephaniehicks/benchmark-hdf5-clustering/blob/2020-05-07-freeze/scripts/simulate_gauss_mix_k.R
simulate_gauss_mix_k <- function(n_cells, n_genes,
                                 k, x_mus, x_sds, 
                                 y_mus, y_sds, prop1)
{ 
  
  if(k != length(x_mus)){stop("k is not same as length of x_mus")} 
  if(k != length(x_sds)){stop("k is not same as length of x_sds")} 
  if(k != length(y_mus)){stop("k is not same as length of y_mus")} 
  if(k != length(y_sds)){stop("k is not same as length of y_sds")} 
  if(k != length(prop1)){stop("k is not same as length of prop1")} 
  
  comp1 <- sample(seq_len(k), prob=prop1, size=n_cells, replace=TRUE)
  
  # Sampling locations for cells in each component
  samples1 <- cbind(rnorm(n=n_cells, mean=x_mus[comp1],sd=x_sds[comp1]),
                    rnorm(n=n_cells, mean=y_mus[comp1],sd=y_sds[comp1]))
  
  # Random projection to D dimensional space, to mimic high-dimensional expression data.
  proj <- matrix(rnorm(n_genes*n_cells), nrow=n_genes, ncol=2)
  A1 <- samples1 %*% t(proj)
  
  # Add normally distributed noise.
  A1 <- A1 + rnorm(n_genes*n_cells)
  rownames(A1) <- paste0("Cell", seq_len(n_cells), "-1")
  colnames(A1) <- paste0("Gene", seq_len(n_genes))
  
  list("true_center" = cbind("x" = x_mus, "y" = y_mus),
       "true_cluster_id" = comp1,
       "true_data" = samples1, 
       "obs_data" = A1)
}

simulate_data <- function(balance, spread, n_cells = 1000, n_genes = 3,
                          k = 3){ 
  if(spread =='easy'){
    x_mus = c(1 ,5, 10)
    y_mus = c(1,5,10)
    x_sds = c(0.1,0.1,0.1)
    y_sds = c(0.1,0.1,0.1)
  } else{ 
    x_mus = c(1 ,2, 3)
    y_mus = c(1,2,3)
    x_sds = c(0.5,0.5,0.5)
    y_sds = c(0.5,0.5,0.5)
  }
  
  if(balance == 'balanced'){
    prop1 = c(0.3, 0.4, 0.3)
  } else{
    prop1 = c(0.1, 0.1, 0.8)
  }
  
  simulate_gauss_mix_k(n_cells = n_cells, n_genes = n_genes, k=k,
                       x_mus = x_mus, y_mus = y_mus, x_sds = x_sds,
                       y_sds = y_sds, prop1 = prop1
  )
}
