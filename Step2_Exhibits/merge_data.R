rm(list=ls())
library(dplyr)
library(ggplot2)
library(R.matlab)
library(zoo)
library(reshape2)
library(abind)
library(RColorBrewer)
library(moments)
library(abind)
library(sandwich)
library(lmtest)
library(openxlsx)

demeanX <- 0
demeanY <- 0

winlist <-c(12,60,120)
for (win in winlist){
  setwd(paste0("~Step1_Predictions/rff_SeparateSims/maxP-12000-trnwin-",win,"-gamma-2-stdize-1-demeanX-",demeanX,"-demeanY-",demeanY,"/"))
  
  # fixed lambda
  ## load information
  mat_files <- list.files(pattern = "^iSim\\d+\\.mat$")
  info <- readMat(mat_files[1]) 
  Y <- as.vector(info$Y)
  D <- as.vector(info$Plist)
  P <- as.vector(info$lamlist)
  dates <- as.vector(info$dates)
  year <- as.numeric(substr(dates, 1, 4))
  n <- length(Y)
  
  ## history strategy
  Yprd_history <- lag(rollapply(Y, width = win, FUN = mean, align = "right", fill = NA))
  timing_history <- Yprd_history*Y
  timing_history <- (1-0.0025*abs(Yprd_history-lag(Yprd_history)))*(1+timing_history) - 1
  
  ## load simple model
  GW <- readMat(paste0("../gybench-trnwin-",win,"-stdize-1-demeanX-",demeanX,"-demeanY-",demeanY,".mat"))
  Yprd_GW <- GW$Yprd.gy
  timing_GW <- Yprd_GW*Y
  timing_GW <- (1-0.0025*abs(Yprd_GW-lag(Yprd_GW)))*(1+timing_GW) - 1
  Hat_GW <- GW$Hat.gy
  
  ## load complex model
  ### read prediction from each sim
  list_Yprd <- list()
  list_Hat <- list()
  file_index <- 0
  for (i in seq_along(mat_files)) {
    file_index <- file_index + 1
    list_Yprd[[i]] <- readMat(mat_files[i])$Yprd
    list_Hat[[i]] <- readMat(mat_files[i])$Hat
    cat(win," train window ",file_index,"th simulation loading...\n")
  }
  ### mean predictions
  Yprd_mean <- Reduce("+", lapply(list_Yprd, function(m) replace(m, is.na(m), 0))) / 
    Reduce("+", lapply(list_Yprd, function(m) !is.na(m)))
  timing_mean <- sweep(Yprd_mean, 1, Y, "*")
  timing_mean <- (1-0.0025*abs(Yprd_mean-lag(Yprd_mean)))*(1+timing_mean) - 1
  Hat_mean <- Reduce("+", list_Hat) / length(list_Hat)
  
  # optimal lambda
  ## load information
  setwd(paste0("~Step1_Predictions/rff_SeparateSims/best-maxP-12000-trnwin-",win,"-gamma-2-stdize-1-demeanX-",demeanX,"-demeanY-",demeanY,"/"))
  mat_files <- list.files(pattern = "^iSim\\d+\\.mat$")
  info <- readMat(mat_files[1]) 
  Y <- as.vector(info$Y)
  D <- as.vector(info$Plist)
  P <- as.vector(info$lamlist)
  dates <- as.vector(info$dates)
  year <- as.numeric(substr(dates, 1, 4))
  n <- length(Y)
  
  ## load optimal lambda for simple model
  all_GW <- readMat(paste0("../best-gybench-trnwin-",win,"-stdize-1-demeanX-",demeanX,"-demeanY-",demeanY,".mat"))
  all_Yprd_GW <- all_GW$Yprd.gy
  all_Hat_GW <- all_GW$Hat
  ### get best lambda for simple model
  best_lambda_GW <- c()
  best_Yprd_GW <- c()
  best_Hat_GW <- c()
  for (t in c((win+1):n)){
    if (t == win + 1){
      best_lambda_idx <- 88 # set lambda to be 10^-3 for first prediction
    } else{
      errors <- apply(all_Yprd_GW[(t-12):(t-1),] - Y[(t-12):(t-1)],2, function(x) sum(x^2,na.rm = TRUE))
      best_lambda_idx <- which.min(errors) # using the minimum error in the lase prediction
    }
    best_lambda_GW[t] <- P[best_lambda_idx]
    best_Yprd_GW[t] <- all_Yprd_GW[t, best_lambda_idx]
    best_Hat_GW[t] <- all_Hat_GW[t, best_lambda_idx]
  }
  best_timing_GW <- best_Yprd_GW*Y
  best_timing_GW <- (1-0.0025*abs(best_Yprd_GW-lag(best_Yprd_GW)))*(1+best_timing_GW) - 1
  
  ## load best lambda for complex model
  list_best_lambda <- list()
  list_best_Yprd <- list()
  list_best_Hat <- list()
  file_index <- 0
  for (i in seq_along(mat_files)) {
    file_index <- file_index + 1
    Yprd <- readMat(mat_files[i])$Yprd
    Hat <- readMat(mat_files[i])$Hat
    list_best_lambda[[i]] = array(NA,c(n,length(D)))
    list_best_Yprd[[i]] = array(NA,c(n,length(D)))
    list_best_Hat[[i]] = array(NA,c(n,length(D)))
    for (t in c((win+1):n)){
      if (t <= win+12){
        best_lambda_idx <- rep(88,length(D))
      } else{
        errors <- apply(Yprd[(t-12):(t-1), , ] - Y[(t-12):(t-1)], c(2,3), function(x) sum(x^2,na.rm = TRUE)) 
        best_lambda_idx <- apply(errors, 1, function(x) {if (all(is.na(x))) 88 else which.min(x)})
      }
      list_best_lambda[[i]][t,] <- P[best_lambda_idx]
      list_best_Yprd[[i]][t,] <- Yprd[t,,][cbind(1:length(D), best_lambda_idx)]
      list_best_Hat[[i]][t,] <- Hat[t,,][cbind(1:length(D), best_lambda_idx)]
    }
    cat(win," train window ",file_index,"th simulation (best) loading...\n")
  }
  best_Yprd_mean <- Reduce("+", lapply(list_best_Yprd, function(m) replace(m, is.na(m), 0))) / 
    Reduce("+", lapply(list_best_Yprd, function(m) !is.na(m)))
  best_timing_mean <- sweep(best_Yprd_mean, 1, Y, "*")
  best_timing_mean <- (1-0.0025*abs(best_Yprd_mean-lag(best_Yprd_mean)))*(1+best_timing_mean) - 1
  best_Hat_mean <- Reduce("+", list_best_Hat) / length(list_best_Hat)
  
  # save everything to .rdata
  save(Y,year,D,
       Yprd_history,timing_history,
       Yprd_GW,best_Yprd_GW,timing_GW,best_timing_GW,Hat_GW,best_Hat_GW,
       Yprd_mean,best_Yprd_mean,timing_mean,best_timing_mean,Hat_mean,best_Hat_mean,
       file = paste0("../../../Step2_Exhibits/",win,"_demeanX-",demeanX,"-demeanY-",demeanY,".RData"))
}