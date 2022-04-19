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

## output directories
dir_plots <- here::here("simulation_figures", "01_across_r" )
dir.create(dir_plots, showWarnings = FALSE)

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

#easy data, balanced 

balance <- c('balanced', 'imbalanced')
spread_levels <- c('easy', 'hard')


easy_data <- simulate_gauss_mix_k(n_cells = 1000, n_genes = 3,
                                  k =3, x_mus = c(1 ,5, 10),
                                  y_mus = c(1,5,10),
                                  x_sds = c(0.1,0.1,0.1),
                                  y_sds = c(0.1,0.1,0.1),
                                  prop1 = c(0.3, 0.4, 0.3))

ks<- seq(4)
kmeans_for_ks <- ks |> map(function(x) kmeans(easy_data$obs_data, centers =x))

fhp_res <- kmeans_for_ks |> map(function(x) fasthplus::hpb(D = easy_data$obs_data, 
                                                           L = x$cluster, t = 10, r = 10)) 

p1 <- plot(ks, fhp_res)


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

easy_balanced_data <- simulate_data(balance = 'balanced', spread = 'easy', k = 3)

ks<- 2:6
kmeans_for_ks <- ks |> map(function(x) kmeans(easy_balanced_data$obs_data, centers =x))

#fhp_res <- kmeans_for_ks |> map(function(x) fasthplus::hpb(D = easy_balanced_data$obs_data, 
                                                       #    L = x$cluster, t = 10, r = 10)) 
fhp_res <-kmeans_for_ks |> map(function(x) fasthplus::hpe(D = dist(data$obs_data), 
                                                          L = x$cluster, 
                                                          alg = "grid_search",
                                                          p= 501))
p2 <- plot(ks, fhp_res)





hard_balanced_data <- simulate_data(balance = 'balanced', spread = 'hard', k = 3)
easy_imbalanced_data <- simulate_data(balance = 'imbalanced', spread = 'easy', k = 3)
hard_imbalanced_data <- simulate_data(balance = 'imbalanced', spread = 'hard', k =3)

par(mar=c(1, 1, 1, 1))
#check true data 
plot(easy_balanced_data$true_data[,1], easy_balanced_data$true_data[,2])
plot(hard_imbalanced_data$true_data[,1], hard_imbalanced_data$true_data[,2],
     col = c(1,2,3)[hard_imbalanced_data$true_cluster_id])




plot_fhp_across_k <- function(data, k, t, r){
  ks<- 2:6
  kmeans_for_ks <- ks |> map(function(x) kmeans(data$obs_data, centers =x))
  
  #fhp_res <- kmeans_for_ks |> map(function(x) fasthplus::hpb(D = data$obs_data, 
                                                             #L = x$cluster, t = t, r = r)) 
  fhp_res <-kmeans_for_ks |> map(function(x) fasthplus::hpe(D = dist(data$obs_data), 
                                                            L = x$cluster, 
                                                            alg = "grid_search",
                                                            p= 501))
                                 
  #plot(ks, fhp_res)
  results <- do.call(rbind, Map(data.frame, ks=ks, fhp_res=fhp_res))
  
  #ggplot( data = results, aes(x=ks, y=fhp_res))+
    #geom_line()+
    #geom_point()

}




##plot x-axis different k, y is a single r 
plot_fhp_across_k(easy_balanced_data, k = 6, t = 10, r = 10)
##replicate calculating hpe for a fixed data
results_list <- replicate(100, plot_fhp_across_k(easy_balanced_data, 
                                                  k = 6, t = 10, r = 10))
melt(results_list[1,], id = c("ks", "fhp_res")) |> mutate(ks = as.character(ks)) |>
  group_by(ks) |> summarise(mean = mean(fhp_res))

# ks      mean
# <chr>  <dbl>
#   2     0.104 
#   3     0.0214
#   4     0.0444
#   5     0.0529
#   6     0.0615


# fit <- kmeans(hard_balanced_data$obs_data, centers = 3)
# fasthplus::hpe(D = dist(hard_balanced_data$obs_data), 
#                L = fit$cluster, p = 101)
# 
# pc <-prcomp(easy_balanced_data$obs_data)


plot_fhp_across_k(easy_imbalanced_data, k = 6, t = 10, r = 10)
plot_fhp_across_k(hard_balanced_data, k = 6, t = 10, r = 10)
plot_fhp_across_k(hard_imbalanced_data, k = 6, t = 10, r = 10)



##plot x- axis different rs, y compare to true H+ 
plot_fhp_across_r<- function(data, k, t, r){
  rs<- seq(r)
  k_means <- kmeans(data$obs_data, centers =k)
  
  fhp_res <- rs |> map(function(x) fasthplus::hpb(D = data$obs_data, 
                                                             L = k_means$cluster,
                                                  t = t, r = x)) 
  
  results <- do.call(rbind, Map(data.frame, rs=rs, fhp_res=fhp_res))
  
  ##compute hpe value 
  hpe_res <- fasthplus::hpe(D=dist(data$obs_data), L = k_means$cluster, p = 101)
  
  ggplot( data = results, aes(x=rs, y=fhp_res))+
    geom_line()+
    geom_point()+
    ##add hpe_res value 
    geom_hline(yintercept = hpe_res)
  
}

plot_fhp_across_r(easy_balanced_data, k = 6, t = 10, r = 10)
plot_fhp_across_r(easy_imbalanced_data, k = 6, t = 10, r = 10)
plot_fhp_across_r(hard_balanced_data, k = 6, t = 10, r = 10)
plot_fhp_across_r(hard_imbalanced_data, k = 6, t = 10, r = 10)


##create plot pdfs 
# all_data<- c(easy_balanced_data,easy_imbalanced_data,hard_balanced_data,
#              hard_imbalanced_data)
# dir_plots <-
#   here::here("plots", "01_across_k", type)
# dir.create(dir_plots, showWarnings = FALSE)
# 



m<-100
n <- 500
cl1 <- sapply(1:n, function(i) rnorm(n=m,mean=0.5,sd=1))
cl2 <- sapply(1:n, function(i) rnorm(n=m,mean=-0.5,sd=1))
cl3 <- sapply(1:n, function(i) rnorm(n=m,mean=-0.5,sd=1))
dat <- t(cbind(cl1,cl2,cl3))
d <- dist(dat)
l <- c(rep(0,n),rep(1,n))
pc <- prcomp(dat)$x[,1:3]
cols <- ifelse(l==1,'#0000ff64','#ff000064')
plot(x=pc[,1],y=pc[,2],pch=16,col=cols,cex=0.7,xaxs = "i",yaxs = "i",xlab='PC1',ylab='PC2',xaxt='n',yaxt='n')
#legend('top',legend=c(expression('Cl'[1]), expression('Cl'[2])), pch=c(21,21),
       #col= c('blue','red'),cex=1.5, pt.bg=c('#0000ff64','#ff000064'),bty='n')

ks<- 2:6
kmeans_for_ks <- ks |> map(function(x) kmeans(dat, centers =x))
fhp_res <- kmeans_for_ks |> map(function(x) fasthplus::hpe(D = dist(dat), 
                                                           L = x$cluster,
                                                           alg = "grid_search",
                                                           p= 1001)) 
fasthplus::hpe(D=d,L=l,p=1001) # D,L formulation

##TODO 
##debug simulated data, explore hpe in easy and hard settings 
