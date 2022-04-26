#!/usr/bin/env Rscript

library(fasthplus)
library(here)
library(usethis)
library(purrr)
library(ggplot2)
library(dplyr)
library(reshape2)

#call external scripts
source(file.path(here::here("simulations","bin","01_simulate_data.R")))

## output directories
dir_plots <- here::here("simulation_figures", "02_hpe_hpb_diff" )
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)


easy_balanced_data <- simulate_data(balance = 'balanced', spread = "easy")


par(mar=c(1, 1, 1, 1))


hpb_across_k <- function(data,  k, t, r){
  ks<- 2:6
  kmeans_for_ks <- ks |> map(function(x) kmeans(data$obs_data, centers =x))
  
  hpb_res <- kmeans_for_ks |> map(function(x) fasthplus::hpb(D = data$obs_data, 
                                                               L = x$cluster, t = t, r = r)) 
  do.call(rbind, Map(data.frame, ks=ks, hpb_res=hpb_res))
  
  }

hpb_across_k(data = easy_balanced_data, t= 100, r = 30)
  
hpe_across_k <- function(data, alg = "grid_search", p = 501){ 
  ks<- 2:6
  kmeans_for_ks <- ks |> map(function(x) kmeans(data$obs_data, centers =x))
  
  hpe_res <-kmeans_for_ks |> map(function(x) fasthplus::hpe(D = dist(data$obs_data), 
                                                            L = x$cluster, 
                                                            alg = alg,
                                                            p= p))
  do.call(rbind, Map(data.frame, ks=ks, hpe_res=hpe_res))
  }
  
  
  
  #ggplot( data = results, aes(x=ks, y=fhp_res))+
    #geom_line()+
    #geom_point()


  

# 
# ##plot x-axis different k, y is a single r 
# plot_fhp_across_k(easy_balanced_data, k = 6, t = 10, r = 10)
# ##replicate calculating hpe for a fixed data
# results_list <- replicate(100, plot_fhp_across_k(easy_balanced_data, 
#                                                   k = 6, t = 10, r = 10))
# melt(results_list[1,], id = c("ks", "fhp_res")) |> mutate(ks = as.character(ks)) |>
#   group_by(ks) |> summarise(mean = mean(fhp_res))

# ks      mean
# <chr>  <dbl>
#   2     0.104 
#   3     0.0214
#   4     0.0444
#   5     0.0529
#   6     0.0615




##plot x- axis different rs, y compare to true H+ 

fhp_across_r<- function(data, k, t = 100, rs = c(10, 50, 100, 200)){
  
  k_means <- kmeans(data$obs_data, centers =k)
  
  fhp_res <- rs |> map(function(x) fasthplus::hpb(D = data$obs_data, 
                                                             L = k_means$cluster,
                                                  t = t, r = x)) 
  
  do.call(rbind, Map(data.frame, rs=rs, fhp_res=fhp_res))
  
}

res <- list()
for(k in ks){
  fhp_across_r()
}



# hpe_res <- fasthplus::hpe(D=dist(data$obs_data), L = k_means$cluster, p = 101)
# 
# ggplot( data = results, aes(x=rs, y=fhp_res))+
#   geom_line()+
#   geom_point()+
#   ##add hpe_res value 
#   geom_hline(yintercept = hpe_res)


