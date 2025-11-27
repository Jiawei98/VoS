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
library(R.matlab)
colors <- brewer.pal(n = 10, name = "Paired")
column_names <- c("beg", "end", "method", "lambda",
                  "R2_1", "SR_1", "MSE_1", "cor_hist_1", "fit1_alpha_t_1","fit2_alpha_t_1", "fit3_alpha_t_1", "fit3_alpha_p_1",
                  "R2_2", "SR_2", "MSE_2", "cor_hist_2", "fit1_alpha_t_2","fit2_alpha_t_2", "fit3_alpha_t_2", "fit3_alpha_p_2")
# 1 for demean, 2 for nodemean
lag = 4
prewhite = FALSE
lam_list = c(-4:4, 9, "best")
linear_ctrn = "3"
setwd("~/Step2_Exhibits/")

winlist <-c(12,60,120)
for (win in winlist){
  table1 <- data.frame(matrix(ncol = length(column_names), nrow = 0))
  table2 <- data.frame(matrix(ncol = length(column_names), nrow = 0))
  
  # get the data
  load(paste0(win,"_demeanX-0-demeanY-0.RData"))
  #simple
  Hat_GW_nodemean <- cbind(Hat_GW,best_Hat_GW)
  colnames(Hat_GW_nodemean) <- lam_list
  timing_GW_nodemean <- cbind(timing_GW,best_timing_GW)
  colnames(timing_GW_nodemean) <- lam_list
  Yprd_GW_nodemean <- cbind(Yprd_GW,best_Yprd_GW)
  colnames(Yprd_GW_nodemean) <- lam_list
  #complex
  Hat_mean_nodemean <- cbind(Hat_mean[,dim(timing_mean)[2],],best_Hat_mean[,dim(timing_mean)[2]])
  colnames(Hat_mean_nodemean) <- lam_list
  timing_mean_nodemean <- cbind(timing_mean[,dim(timing_mean)[2],],best_timing_mean[,dim(timing_mean)[2]])
  colnames(timing_mean_nodemean) <- lam_list
  Yprd_mean_nodemean <- cbind(Yprd_mean[,dim(timing_mean)[2],],best_Yprd_mean[,dim(timing_mean)[2]])
  colnames(Yprd_mean_nodemean) <- lam_list
  load(paste0(win,"_demeanX-1-demeanY-1.RData"))
  #simple
  Hat_GW_demean <- cbind(Hat_GW,best_Hat_GW)
  colnames(Hat_GW_demean) <- lam_list
  timing_GW_demean <- cbind(timing_GW,best_timing_GW)
  colnames(timing_GW_demean) <- lam_list
  Yprd_GW_demean <- cbind(Yprd_GW,best_Yprd_GW)
  colnames(Yprd_GW_demean) <- lam_list
  #complex
  Hat_mean_demean <- cbind(Hat_mean[,dim(timing_mean)[2],],best_Hat_mean[,dim(timing_mean)[2]])
  colnames(Hat_mean_demean) <- lam_list
  timing_mean_demean <- cbind(timing_mean[,dim(timing_mean)[2],],best_timing_mean[,dim(timing_mean)[2]])
  colnames(timing_mean_demean) <- lam_list
  Yprd_mean_demean <- cbind(Yprd_mean[,dim(timing_mean)[2],],best_Yprd_mean[,dim(timing_mean)[2]])
  colnames(Yprd_mean_demean) <- lam_list
  
  # split data
  subbeg <- c(1930,1975,1930)
  subend <- c(1974,2020,2020)
  
  # subsample
  for (t in 1:length(subbeg)){
    indices <- !is.na(timing_history) & year >= subbeg[t] & year <= subend[t]
    n0 <- sum(indices)
    
    ### get the predictions
    Y_temp <- Y[indices]
    timing_history_temp <- timing_history[indices]
    Yprd_history_temp <- Yprd_history[indices]
    timing_GW_nodemean_temp <- timing_GW_nodemean[indices,]
    Yprd_GW_nodemean_temp <- Yprd_GW_nodemean[indices,]
    timing_mean_nodemean_temp <- timing_mean_nodemean[indices,]
    Yprd_mean_nodemean_temp <- Yprd_mean_nodemean[indices,]
    timing_GW_demean_temp <- timing_GW_demean[indices,]
    Yprd_GW_demean_temp <- Yprd_GW_demean[indices,]
    timing_mean_demean_temp <- timing_mean_demean[indices,]
    Yprd_mean_demean_temp <- Yprd_mean_demean[indices,]

    ### MKT
    SR <- sqrt(12) * mean(Y_temp, na.rm = TRUE) / sd(Y_temp, na.rm = TRUE) # * sqrt(n0/(n0-1))
    
    ##### add to table
    table1 <- rbind(table1,
                    list(subbeg[t], subend[t], "market", NA,
                         NA, SR, NA, NA, NA, NA, NA, NA,
                         NA, SR, NA, NA, NA, NA, NA, NA))
    table2 <- rbind(table2,
                    list(subbeg[t], subend[t], "market", NA,
                         NA, SR, NA, NA, NA, NA, NA, NA,
                         NA, SR, NA, NA, NA, NA, NA, NA))
    
    ### history
    R2 <- 1 - sum((Yprd_history_temp - Y_temp)^2, na.rm = TRUE) / sum((Y_temp - Yprd_history_temp)^2, na.rm = TRUE)
    MSE <- sum((Yprd_history_temp - Y_temp)^2, na.rm = TRUE)
    SR <- sqrt(12) * mean(timing_history_temp, na.rm = TRUE) / sd(timing_history_temp, na.rm = TRUE) # * sqrt(n0/(n0-1))
    
    ##### fit1 timing_history ~ MKT
    fit <- lm(timing_history_temp ~ Y_temp)
    NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
    stats <- coeftest(fit, vcov = NW_VCOV)
    fit1_alpha <- stats[1, "Estimate"]
    fit1_alpha_t <- stats[1, "t value"]
    fit1_alpha_p <- stats[1, "Pr(>|t|)"]
    
    ##### add to table
    table1 <- rbind(table1,
                    list(subbeg[t], subend[t], "history", NA,
                         R2, SR, MSE, NA, fit1_alpha_t, NA, NA, NA,
                         R2, SR, MSE, NA, fit1_alpha_t, NA, NA, NA))
    table2 <- rbind(table2,
                    list(subbeg[t], subend[t], "history", NA,
                         R2, SR, MSE, NA, fit1_alpha_t, NA, NA, NA,
                         R2, SR, MSE, NA, fit1_alpha_t, NA, NA, NA))
    ### simple model
    for (lambda in lam_list){
      #### demean
      ##### performance metrics
      demean_R2 <- 1 - sum((Yprd_GW_demean_temp[,lambda] - Y_temp)^2, na.rm = TRUE) / sum((Y_temp - Yprd_history_temp)^2, na.rm = TRUE)
      demean_MSE <- sum((Yprd_GW_demean_temp[,lambda] - Y_temp)^2, na.rm = TRUE)
      demean_SR <- sqrt(12) * mean(timing_GW_demean_temp[,lambda], na.rm = TRUE) / sd(timing_GW_demean_temp[,lambda], na.rm = TRUE) # * sqrt(n0/(n0-1))
      demean_cor <- cor(timing_GW_demean_temp[,lambda],timing_history_temp)
      
      ##### fit1 timing_GW ~ MKT
      fit <- lm(timing_GW_demean_temp[,lambda] ~ Y_temp)
      NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
      stats <- coeftest(fit, vcov = NW_VCOV)
      if (lag == 0) {stats <- summary(fit)$coefficients}
      demean_fit1_alpha <- stats[1, "Estimate"]
      demean_fit1_alpha_t <- stats[1, "t value"]
      demean_fit1_alpha_p <- stats[1, "Pr(>|t|)"]
      
      ##### fit2 timing_GW ~ MKT + timing_history
      fit <- lm(timing_GW_demean_temp[,lambda] ~ Y_temp+timing_history_temp)
      NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
      stats <- coeftest(fit, vcov = NW_VCOV)
      if (lag == 0) {stats <- summary(fit)$coefficients}
      demean_fit2_alpha <- stats[1, "Estimate"]
      demean_fit2_alpha_t <- stats[1, "t value"]
      demean_fit2_alpha_p <- stats[1, "Pr(>|t|)"]
      
      #### nodemean
      ##### performance metrics
      nodemean_R2 <- 1 - sum((Yprd_GW_nodemean_temp[,lambda] - Y_temp)^2, na.rm = TRUE) / sum((Y_temp - Yprd_history_temp)^2, na.rm = TRUE)
      nodemean_MSE <- sum((Yprd_GW_nodemean_temp[,lambda] - Y_temp)^2, na.rm = TRUE)
      nodemean_SR <- sqrt(12) * mean(timing_GW_nodemean_temp[,lambda], na.rm = TRUE) / sd(timing_GW_nodemean_temp[,lambda], na.rm = TRUE) # * sqrt(n0/(n0-1))
      nodemean_cor <- cor(timing_GW_nodemean_temp[,lambda],timing_history_temp)
      
      ##### fit1 timing_GW ~ MKT
      fit <- lm(timing_GW_nodemean_temp[,lambda] ~ Y_temp)
      NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
      stats <- coeftest(fit, vcov = NW_VCOV)
      if (lag == 0) {stats <- summary(fit)$coefficients}
      nodemean_fit1_alpha <- stats[1, "Estimate"]
      nodemean_fit1_alpha_t <- stats[1, "t value"]
      nodemean_fit1_alpha_p <- stats[1, "Pr(>|t|)"]
      
      ##### fit2 timing_GW ~ MKT + timing_history
      fit <- lm(timing_GW_nodemean_temp[,lambda] ~ Y_temp+timing_history_temp)
      NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
      stats <- coeftest(fit, vcov = NW_VCOV)
      if (lag == 0) {stats <- summary(fit)$coefficients}
      nodemean_fit2_alpha <- stats[1, "Estimate"]
      nodemean_fit2_alpha_t <- stats[1, "t value"]
      nodemean_fit2_alpha_p <- stats[1, "Pr(>|t|)"]
      
      ##### add to table
      table1 <- rbind(table1,
                      list(subbeg[t], subend[t], "linear", lambda, 
                           demean_R2,demean_SR,demean_MSE,demean_cor,demean_fit1_alpha_t,demean_fit2_alpha_t, NA, NA,
                           nodemean_R2,nodemean_SR,nodemean_MSE,nodemean_cor,nodemean_fit1_alpha_t,nodemean_fit2_alpha_t, NA, NA))
      if (lambda == "3"){
        table2 <- rbind(table2,
                        list(subbeg[t], subend[t], "linear", lambda, 
                             demean_R2,demean_SR,demean_MSE,demean_cor,demean_fit1_alpha_t,demean_fit2_alpha_t, NA, NA,
                             nodemean_R2,nodemean_SR,nodemean_MSE,nodemean_cor,nodemean_fit1_alpha_t,nodemean_fit2_alpha_t, NA, NA))
      }
    }
    
    ### complex model
    for (lambda in lam_list){
      #### demean
      ##### performance metrics
      demean_R2 <- 1 - sum((Yprd_mean_demean_temp[,lambda] - Y_temp)^2, na.rm = TRUE) / sum((Y_temp - Yprd_history_temp)^2, na.rm = TRUE)
      demean_MSE <- sum((Yprd_mean_demean_temp[,lambda] - Y_temp)^2, na.rm = TRUE)
      demean_SR <- sqrt(12) * mean(timing_mean_demean_temp[,lambda], na.rm = TRUE) / sd(timing_mean_demean_temp[,lambda], na.rm = TRUE) # * sqrt(n0/(n0-1))
      demean_cor <- cor(timing_mean_demean_temp[,lambda],timing_history_temp)
      
      ##### fit1 timing_GW ~ MKT
      fit <- lm(timing_mean_demean_temp[,lambda] ~ Y_temp)
      NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
      stats <- coeftest(fit, vcov = NW_VCOV)
      if (lag == 0) {stats <- summary(fit)$coefficients}
      demean_fit1_alpha <- stats[1, "Estimate"]
      demean_fit1_alpha_t <- stats[1, "t value"]
      demean_fit1_alpha_p <- stats[1, "Pr(>|t|)"]
      
      ##### fit2 timing_GW ~ MKT + timing_history
      fit <- lm(timing_mean_demean_temp[,lambda] ~ Y_temp+timing_history_temp)
      NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
      stats <- coeftest(fit, vcov = NW_VCOV)
      if (lag == 0) {stats <- summary(fit)$coefficients}
      demean_fit2_alpha <- stats[1, "Estimate"]
      demean_fit2_alpha_t <- stats[1, "t value"]
      demean_fit2_alpha_p <- stats[1, "Pr(>|t|)"]
      
      ##### fit3 timing_GW ~ MKT + timing_history + demean_simple
      if (t == 1 | t == 2){
        fit <- lm(timing_mean_demean_temp[,lambda] ~ Y_temp+timing_history_temp+timing_GW_demean_temp[,linear_ctrn])
        NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
        stats <- coeftest(fit, vcov = NW_VCOV)
        if (lag == 0) {stats <- summary(fit)$coefficients}
        demean_fit3_alpha <- stats[1, "Estimate"]
        demean_fit3_alpha_t <- stats[1, "t value"]
        demean_fit3_alpha_p <- stats[1, "Pr(>|t|)"]
      } else {
        # fit first half
        indices <- !is.na(Yprd_GW[,1]) & year >= subbeg[1] & year <= subend[1]
        fit1 <- lm(timing_mean_demean[indices,lambda] ~ Y[indices]+timing_history[indices]+timing_GW_demean[indices,linear_ctrn])
        NW_VCOV <- NeweyWest(fit1, lag = lag, prewhite = prewhite)
        stats <- coeftest(fit1, vcov = NW_VCOV)
        if (lag == 0) {stats <- summary(fit)$coefficients}
        alpha <- stats[1, "Estimate"]
        res1 <- alpha + residuals(fit1)
        # fit second half
        indices <- !is.na(Yprd_GW[,1]) & year >= subbeg[2] & year <= subend[2]
        fit2 <- lm(timing_mean_demean[indices,lambda] ~ Y[indices]+timing_history[indices]+timing_GW_demean[indices,linear_ctrn])
        NW_VCOV <- NeweyWest(fit2, lag = lag, prewhite = prewhite)
        stats <- coeftest(fit2, vcov = NW_VCOV)
        if (lag == 0) {stats <- summary(fit)$coefficients}
        alpha <- stats[1, "Estimate"]
        res2 <- alpha + residuals(fit2)
        # fit full
        res <- c(res1, res2)
        fit <- lm(res ~ 1)
        NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
        stats <- coeftest(fit, vcov = NW_VCOV)
        if (lag == 0) {stats <- summary(fit)$coefficients}
        demean_fit3_alpha <- stats[1, "Estimate"]
        demean_fit3_alpha_t <- stats[1, "t value"]
        demean_fit3_alpha_p <- stats[1, "Pr(>|t|)"]
      }
      
      #### nodemean
      ##### performance metrics
      nodemean_R2 <- 1 - sum((Yprd_mean_nodemean_temp[,lambda] - Y_temp)^2, na.rm = TRUE) / sum((Y_temp - Yprd_history_temp)^2, na.rm = TRUE)
      nodemean_MSE <- sum((Yprd_mean_nodemean_temp[,lambda] - Y_temp)^2, na.rm = TRUE)
      nodemean_SR <- sqrt(12) * mean(timing_mean_nodemean_temp[,lambda], na.rm = TRUE) / sd(timing_mean_nodemean_temp[,lambda], na.rm = TRUE) # * sqrt(n0/(n0-1))
      nodemean_cor <- cor(timing_mean_nodemean_temp[,lambda],timing_history_temp)
      
      ##### fit1 timing_GW ~ MKT
      fit <- lm(timing_mean_nodemean_temp[,lambda] ~ Y_temp)
      NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
      stats <- coeftest(fit, vcov = NW_VCOV)
      if (lag == 0) {stats <- summary(fit)$coefficients}
      nodemean_fit1_alpha <- stats[1, "Estimate"]
      nodemean_fit1_alpha_t <- stats[1, "t value"]
      nodemean_fit1_alpha_p <- stats[1, "Pr(>|t|)"]
      
      ##### fit2 timing_GW ~ MKT + timing_history
      fit <- lm(timing_mean_nodemean_temp[,lambda] ~ Y_temp+timing_history_temp)
      NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
      stats <- coeftest(fit, vcov = NW_VCOV)
      if (lag == 0) {stats <- summary(fit)$coefficients}
      nodemean_fit2_alpha <- stats[1, "Estimate"]
      nodemean_fit2_alpha_t <- stats[1, "t value"]
      nodemean_fit2_alpha_p <- stats[1, "Pr(>|t|)"]
      
      ##### fit3 timing_GW ~ MKT + timing_history + demean_simple
      if (t == 1 | t == 2){
        fit <- lm(timing_mean_nodemean_temp[,lambda] ~ Y_temp+timing_history_temp+timing_GW_demean_temp[,linear_ctrn]+timing_GW_nodemean_temp[,linear_ctrn])
        NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
        stats <- coeftest(fit, vcov = NW_VCOV)
        if (lag == 0) {stats <- summary(fit)$coefficients}
        nodemean_fit3_alpha <- stats[1, "Estimate"]
        nodemean_fit3_alpha_t <- stats[1, "t value"]
        nodemean_fit3_alpha_p <- stats[1, "Pr(>|t|)"]
      } else {
        # fit first half
        indices <- !is.na(Yprd_GW[,1]) & year >= subbeg[1] & year <= subend[1]
        fit1 <- lm(timing_mean_nodemean[indices,lambda] ~ Y[indices]+timing_history[indices]+timing_GW_demean[indices,linear_ctrn]+timing_GW_nodemean[indices,linear_ctrn])
        NW_VCOV <- NeweyWest(fit1, lag = lag, prewhite = prewhite)
        stats <- coeftest(fit1, vcov = NW_VCOV)
        if (lag == 0) {stats <- summary(fit)$coefficients}
        alpha <- stats[1, "Estimate"]
        res1 <- alpha + residuals(fit1)
        # fit second half
        indices <- !is.na(Yprd_GW[,1]) & year >= subbeg[2] & year <= subend[2]
        fit2 <- lm(timing_mean_nodemean[indices,lambda] ~ Y[indices]+timing_history[indices]+timing_GW_demean[indices,linear_ctrn]+timing_GW_nodemean[indices,linear_ctrn])
        NW_VCOV <- NeweyWest(fit2, lag = lag, prewhite = prewhite)
        stats <- coeftest(fit2, vcov = NW_VCOV)
        if (lag == 0) {stats <- summary(fit)$coefficients}
        alpha <- stats[1, "Estimate"]
        res2 <- alpha + residuals(fit2)
        # fit full
        res <- c(res1, res2)
        fit <- lm(res ~ 1)
        NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
        stats <- coeftest(fit, vcov = NW_VCOV)
        if (lag == 0) {stats <- summary(fit)$coefficients}
        nodemean_fit3_alpha <- stats[1, "Estimate"]
        nodemean_fit3_alpha_t <- stats[1, "t value"]
        nodemean_fit3_alpha_p <- stats[1, "Pr(>|t|)"]
      }
      
      ##### add to table
      table1 <- rbind(table1,
                      list(subbeg[t], subend[t], "nonlinear", lambda, 
                           demean_R2,demean_SR,demean_MSE,demean_cor,demean_fit1_alpha_t,demean_fit2_alpha_t,demean_fit3_alpha_t,demean_fit3_alpha_p,
                           nodemean_R2,nodemean_SR,nodemean_MSE,nodemean_cor,nodemean_fit1_alpha_t,nodemean_fit2_alpha_t,nodemean_fit3_alpha_t,nodemean_fit3_alpha_p))
      if (lambda == "3" | lambda == "best"){
        table2 <- rbind(table2,
                        list(subbeg[t], subend[t], "nonlinear", lambda, 
                             demean_R2,demean_SR,demean_MSE,demean_cor,demean_fit1_alpha_t,demean_fit2_alpha_t,demean_fit3_alpha_t,demean_fit3_alpha_p,
                             nodemean_R2,nodemean_SR,nodemean_MSE,nodemean_cor,nodemean_fit1_alpha_t,nodemean_fit2_alpha_t,nodemean_fit3_alpha_t,nodemean_fit3_alpha_p))
      }
    }
    
  }
  colnames(table1) <- colnames(table2) <- column_names
  write.csv(table1,file = paste0(win,"_full.csv"))
  write.csv(table2,file = paste0(win,"_final.csv"))
}