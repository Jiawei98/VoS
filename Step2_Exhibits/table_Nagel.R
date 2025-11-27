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
                  "R2_1", "SR_1", "MSE_1", "cor_hist", 
                  "fit1_alpha_t_1","fit2_alpha_t_1", "fit3_alpha_t_1", "fit3_alpha_p_1")
table1 <- data.frame(matrix(ncol = length(column_names), nrow = 0))

win = 12
lag = 4
prewhite = FALSE
setwd("~/Step2_Exhibits/")
M <- c('Kernel', 'EW', 'DW', 'PV', 'EWPV', 'DWPV', 'EWSD')

# read result
load(paste0(win,"_demeanX-1-demeanY-1.RData"))
Yprd_legal <- readMat("Nagel_12.mat")

# Define the desired new order and names
new_order <- c("yhat.kernel", "yhat.EW", "yhat.DW", "yhat.PV", "yhat.EWPV", "yhat.DWPV", "yhat.EWSD")
new_names <- c("Kernel", "EW", "DW", "PV", "EWPV", "DWPV", "EWSD")

# Reorder and rename
Yprd_legal <- Yprd_legal[new_order]
names(Yprd_legal) <- new_names

# timing stragety
timing_legal <- lapply(Yprd_legal, function(x) x * Y)
timing_legal <- lapply(seq_along(timing_legal), function(i) {
  timing_legali <- c(timing_legal[[i]])
  Yprd_legali <- c(Yprd_legal[[i]])
  (1 - 0.0025 * abs(Yprd_legali - lag(Yprd_legali))) * (1 + timing_legali) - 1
})
timing_GW_103 <- timing_GW[,8]

# split data
subbeg <- c(1930,1975,1930)
subend <- c(1974,2020,2020)

# subsample
for (t in 1:length(subbeg)){
  indices <- !is.na(timing_history) & year >= subbeg[t] & year <= subend[t]
  n0 <- sum(indices)
  
  # filter data
  Y_temp <- Y[indices]
  timing_history_temp <- timing_history[indices]
  timing_GW_103_temp <- timing_GW_103[indices]
  Yprd_history_temp <- Yprd_history[indices]
  Yprd_legal_temp <- lapply(Yprd_legal, function(x) x[indices])
  timing_legal_temp <- lapply(timing_legal, function(x) x[indices])
  
  for (i in 1:length(Yprd_legal_temp)){
    R2 <- 1 - sum((Yprd_legal_temp[[i]] - Y_temp)^2, na.rm = TRUE) / sum((Y_temp - Yprd_history_temp)^2, na.rm = TRUE)
    MSE <- sum((Yprd_legal_temp[[i]] - Y_temp)^2, na.rm = TRUE)
    SR <- sqrt(12) * mean(timing_legal_temp[[i]], na.rm = TRUE) / sd(timing_legal_temp[[i]], na.rm = TRUE)
    cor <- cor(timing_legal_temp[[i]],timing_history_temp)
    
    ##### fit1 timing_GW ~ MKTYes
    fit <- lm(timing_legal_temp[[i]] ~ Y_temp)
    NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
    stats <- coeftest(fit, vcov = NW_VCOV)
    if (lag == 0) {stats <- summary(fit)$coefficients}
    fit1_alpha <- stats[1, "Estimate"]
    fit1_alpha_t <- stats[1, "t value"]
    fit1_alpha_p <- stats[1, "Pr(>|t|)"]
    
    ##### fit2 timing_GW ~ MKT + timing_history
    if (i == 2){
      fit2_alpha <- NA
      fit2_alpha_t <- NA
      fit2_alpha_p <- NA
    } else{
      fit <- lm(timing_legal_temp[[i]] ~ Y_temp+timing_history_temp)
      NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
      stats <- coeftest(fit, vcov = NW_VCOV)
      if (lag == 0) {stats <- summary(fit)$coefficients}
      fit2_alpha <- stats[1, "Estimate"]
      fit2_alpha_t <- stats[1, "t value"]
      fit2_alpha_p <- stats[1, "Pr(>|t|)"]
    }
    
    ##### fit3 timing_GW ~ MKT + timing_history + demean_simple
    if (t == 1 | t == 2){
      if (i == 2){
        fit3_alpha <- NA
        fit3_alpha_t <- NA
        fit3_alpha_p <- NA
      } else{
        fit <- lm(timing_legal_temp[[i]] ~ Y_temp+timing_history_temp+timing_GW_103_temp)
        NW_VCOV <- NeweyWest(fit, lag = lag, prewhite = prewhite)
        stats <- coeftest(fit, vcov = NW_VCOV)
        if (lag == 0) {stats <- summary(fit)$coefficients}
        fit3_alpha <- stats[1, "Estimate"]
        fit3_alpha_t <- stats[1, "t value"]
        fit3_alpha_p <- stats[1, "Pr(>|t|)"]
      }
    } else {
      if (i == 2){
        fit3_alpha <- NA
        fit3_alpha_t <- NA
        fit3_alpha_p <- NA
      } else{
        # fit first half
        indices <- !is.na(Yprd_GW[,1]) & year >= subbeg[1] & year <= subend[1]
        fit1 <- lm(timing_legal[[i]][indices] ~ Y[indices]+timing_history[indices]+timing_GW_103[indices])
        NW_VCOV <- NeweyWest(fit1, lag = lag, prewhite = prewhite)
        stats <- coeftest(fit1, vcov = NW_VCOV)
        if (lag == 0) {stats <- summary(fit)$coefficients}
        alpha <- stats[1, "Estimate"]
        res1 <- alpha + residuals(fit1)
        # fit second half
        indices <- !is.na(Yprd_GW[,1]) & year >= subbeg[2] & year <= subend[2]
        fit2 <- lm(timing_legal[[i]][indices] ~ Y[indices]+timing_history[indices]+timing_GW_103[indices])
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
        fit3_alpha <- stats[1, "Estimate"]
        fit3_alpha_t <- stats[1, "t value"]
        fit3_alpha_p <- stats[1, "Pr(>|t|)"]
      }
    }
    
    ##### add to table
    table1 <- rbind(table1,
                    list(subbeg[t], subend[t], M[i], NA, 
                         R2,SR,MSE,cor,fit1_alpha_t,fit2_alpha_t,fit3_alpha_t,fit3_alpha_p))
  }
}
colnames(table1) <- column_names
write.csv(table1,file = paste0(win,"_Nagel.csv"))
