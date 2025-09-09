#' Format p-value for display
#' @param pval Numeric p-value.
#' @param eps Threshold for small p-values.
#' @param digits Number of digits to display.
#' @return Formatted p-value as character.
#' @export
format_pval <- function(pval, eps = 0.001, digits = 3) {
  if (is.na(pval)) return(NA)
  if (pval < eps) return(paste0("<", eps))
  format(round(pval, digits), nsmall = digits)
}

#' Validate required columns in a data frame
#' @param df Data frame to check.
#' @param required_cols Character vector of required column names.
#' @return NULL if all columns present, otherwise error.
#' @export
validate_input <- function(df, required_cols) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) stop(paste("Missing required columns:", paste(missing, collapse = ", ")))
  invisible(NULL)
}


# Weighted counting process: number of events up to time x
#' Weighted counting process
#' @param x Time point.
#' @param y Vector of event/censoring times.
#' @param w Weights (default 1).
#' @return Weighted count of events up to time x.
#' @export
count_weighted <- function(x, y, w = rep(1, length(y))) {
  sum(w * (y <= x))
}

# Weighted risk set: number at risk at time x
#' Weighted risk set
#' @param x Time point.
#' @param y Vector of event/censoring times.
#' @param w Weights (default 1).
#' @return Weighted number at risk at time x.
#' @export
risk_weighted <- function(x, y, w = rep(1, length(y))) {
  sum(w * (y >= x))
}

#' Kaplan-Meier quantile calculation
#' @param time_points Vector of time points.
#' @param survival_probs Vector of survival probabilities.
#' @param qprob Quantile probability (default 0.5).
#' @param type Calculation type (midpoint or min).
#' @return Estimated quantile time.
#' @export
kmq_calculations <- function(time_points, survival_probs, qprob = 0.5, type = "midpoint") {
  tq2 <- suppressWarnings(min(time_points[which(survival_probs <= qprob)]))
  loc.tq2 <- which(time_points == tq2)
  # jump point prior to quant
  qjt1 <- suppressWarnings(min(survival_probs[survival_probs > qprob]))
  tq1 <- suppressWarnings(min(time_points[which(survival_probs == qjt1)]))
  tq.hat <- tq2
  mid_flag <- !is.na(qjt1) && (round(qjt1, 12) == qprob)
  if (type == "midpoint" && mid_flag) {
    tq.hat <- tq1 + (tq2 - tq1) / 2
  }
  if (is.infinite(tq.hat) || is.na(tq.hat)) tq.hat <- NA
  return(tq.hat)
}

#' Kaplan-Meier quantile and confidence interval
#' @param time_points Vector of time points.
#' @param survival_probs Vector of survival probabilities.
#' @param se_probs Standard errors of survival probabilities.
#' @param qprob Quantile probability (default 0.5).
#' @param type Calculation type (midpoint or min).
#' @param conf_level Confidence level (default 0.95).
#' @return List with quantile and confidence interval.
#' @export
km_quantile <- function(time_points, survival_probs, se_probs = NULL, qprob = 0.5, type = c("midpoint","min"), conf_level = 0.95) {
  type <- match.arg(type)
  z <- qnorm(1 - (1 - conf_level) / 2)
  qhat <- kmq_calculations(time_points = time_points, survival_probs = survival_probs, qprob = qprob, type = type)
  qhat[is.infinite(qhat)] <- NA
  qhat_lower <- qhat_upper <- NA
  if (!is.null(se_probs)) {
    # log transform for CI
    lower_probs <- exp(log(survival_probs) - z * se_probs / survival_probs)
    upper_probs <- exp(log(survival_probs) + z * se_probs / survival_probs)
    qhat_lower <- kmq_calculations(time_points = time_points, survival_probs = lower_probs, qprob = qprob, type = type)
    qhat_lower[is.infinite(qhat_lower)] <- NA
    qhat_upper <- kmq_calculations(time_points = time_points, survival_probs = upper_probs, qprob = qprob, type = type)
    qhat_upper[is.infinite(qhat_upper)] <- NA
  }
  return(list(qhat = qhat, lower = qhat_lower, upper = qhat_upper))
}

#' Table of KM quantiles for two groups
#' @param time_points Vector of time points.
#' @param surv0 Survival probabilities for group 0.
#' @param se0 Standard errors for group 0.
#' @param surv1 Survival probabilities for group 1.
#' @param se1 Standard errors for group 1.
#' @param arms Group labels.
#' @param qprob Quantile probability.
#' @param type Calculation type.
#' @param conf_level Confidence level.
#' @return Data frame of quantiles and CIs for each group.
#' @export
km_quantile_table <- function(time_points, surv0, se0, surv1, se1, arms = c("treat", "control"), qprob = 0.5, type = c("midpoint","min"), conf_level = 0.95) {
  type <- match.arg(type)
  kmq0 <- km_quantile(time_points = time_points, survival_probs = surv0, se_probs = se0, qprob = qprob, type = type, conf_level = conf_level)
  df0 <- data.frame(group = arms[2], quantile = kmq0$qhat, lower = kmq0$lower, upper = kmq0$upper)
  kmq1 <- km_quantile(time_points = time_points, survival_probs = surv1, se_probs = se1, qprob = qprob, type = type, conf_level = conf_level)
  df1 <- data.frame(group = arms[1], quantile = kmq1$qhat, lower = kmq1$lower, upper = kmq1$upper)
  quantiles_df <- rbind(df1, df0)
  return(quantiles_df)
}

#' Kaplan-Meier estimates and Greenwood variance
#' @param ybar Number at risk at each time.
#' @param nbar Number of events at each time.
#' @return List with survival and variance estimates.
#' @export
KM_estimates <- function(ybar, nbar, sig2w_multiplier = NULL){
  dN <- diff(c(0, nbar))
  dN_risk <- ifelse(ybar > 0, dN / ybar, 0.0)
  S_KM <- cumprod(1 - dN_risk)
  if(is.null(sig2w_multiplier)){
  sig2w_multiplier  <- ifelse(ybar > 0 & ybar > dN, dN / (ybar * (ybar - dN)), 0.0)
  }
  var_KM <- (S_KM^2) * cumsum(sig2w_multiplier)
  list(S_KM = S_KM, sig2_KM = var_KM)
}

#' Weighted log-rank estimates and variance
#' @param ybar0 Number at risk in group 0.
#' @param ybar1 Number at risk in group 1.
#' @param nbar0 Number of events in group 0.
#' @param nbar1 Number of events in group 1.
#' @param rho Weighting parameter.
#' @param gamma Weighting parameter.
#' @return List with log-rank statistic and variance.
#' @export
wlr_estimates <- function(ybar0, ybar1, nbar0, nbar1, S_pool = NULL, rho = 0, gamma = 0) {
  dN_z0 <- diff(c(0, nbar0))
  dN_z1 <- diff(c(0, nbar1))
  dN_pooled <- dN_z0 + dN_z1
  risk_z1 <- ybar1
  risk_z0 <- ybar0
  risk_pooled <- risk_z0 + risk_z1
  w <- rep(1, length(ybar0))
  if( rho != 0.0 | gamma != 0.0){
  if(is.null(S_pool)){
  dN_Risk <- ifelse(risk_pooled > 0, dN_pooled / risk_pooled, 0)
  S_pool <- cumprod(1 - dN_Risk)
  S_pool <- c(1, S_pool[-length(S_pool)]) # S_pool(t-)
  }
  if(!is.null(S_pool))  S_pool <- c(1, S_pool[-length(S_pool)])
  w <- (S_pool^rho) * ((1 - S_pool)^gamma)
  }
  K <- ifelse(risk_pooled > 0, w * (risk_z0 * risk_z1) / risk_pooled, 0.0)
  drisk0 <- sum(ifelse(risk_z0 > 0, (K / risk_z0) * dN_z0, 0.0))
  drisk1 <- sum(ifelse(risk_z1 > 0, (K / risk_z1) * dN_z1, 0.0))
  lr <- drisk0 - drisk1
  h0 <- ifelse(risk_z0 == 0, 0, (K^2 / risk_z0))
  h1 <- ifelse(risk_z1 == 0, 0, (K^2 / risk_z1))
  dJ <- ifelse(risk_pooled == 1, 0, (dN_pooled - 1) / (risk_pooled - 1))
  dL <- ifelse(risk_pooled == 0, 0, dN_pooled / risk_pooled)
  sig2 <- sum((h0 + h1) * (1 - dJ) * dL)
  list(lr = lr, sig2 = sig2)
}

#' Weighted log-rank and KM difference at tzero
#' @param dfcounting Data frame with counting process columns.
#' @param rho Weighting parameter.
#' @param gamma Weighting parameter.
#' @param tzero Time point for difference.
#' @return List with statistics, variances, covariance, and correlation.
#' @export
wlr_dhat_estimates <- function(dfcounting, rho = 0, gamma = 0, tzero = 24) {
  at_points <- dfcounting$at.points
  nbar0 <- dfcounting$nbar0
  nbar1 <- dfcounting$nbar1
  ybar0 <- dfcounting$ybar0
  ybar1 <- dfcounting$ybar1
  S1 <- dfcounting$surv1
  S0 <- dfcounting$surv0
  S_pool <- dfcounting$survP
  loc_tzero <- which.max(at_points > tzero)
  if (at_points[loc_tzero] <= tzero & at_points[loc_tzero + 1] > tzero) {
    dhat_tzero <- S1[loc_tzero] - S0[loc_tzero]
  } else {
    dhat_tzero <- S1[loc_tzero - 1] - S0[loc_tzero - 1]
  }
  Sp_tzero <- S_pool[loc_tzero]
  dN_z0 <- diff(c(0, nbar0))
  dN_z1 <- diff(c(0, nbar1))
  dN_pooled <- dN_z0 + dN_z1
  risk_z1 <- ybar1
  risk_z0 <- ybar0
  risk_pooled <- risk_z0 + risk_z1
  S_pool <- c(1, S_pool[-length(S_pool)])
  w <- (S_pool^rho) * ((1 - S_pool)^gamma)
  K <- ifelse(risk_pooled > 0, w * (risk_z0 * risk_z1) / risk_pooled, 0.0)
  drisk0 <- sum(ifelse(risk_z0 > 0, (K / risk_z0) * dN_z0, 0.0))
  drisk1 <- sum(ifelse(risk_z1 > 0, (K / risk_z1) * dN_z1, 0.0))
  lr <- drisk0 - drisk1
  h0 <- ifelse(risk_z0 == 0, 0, (K^2 / risk_z0))
  h1 <- ifelse(risk_z1 == 0, 0, (K^2 / risk_z1))
  dJ <- ifelse(risk_pooled == 1, 0, (dN_pooled - 1) / (risk_pooled - 1))
  dL <- ifelse(risk_pooled == 0, 0, dN_pooled / risk_pooled)
  sig2_lr <- sum((h0 + h1) * (1 - dJ) * dL)
  w_tzero <- w * ifelse(at_points <= tzero, 1, 0)
  w_integral_t0 <- sum(w_tzero * (1 - dJ) * dL)
  cov_wlr_dhat <- Sp_tzero * w_integral_t0
  h2 <- ifelse(risk_z0 * risk_z1 > 0, (risk_pooled / (risk_z0 * risk_z1)), 0)
  h2 <- h2 * ifelse(at_points <= tzero, 1, 0)
  sig2_dhat <- (Sp_tzero^2) * sum(h2 * (1 - dJ) * dL)
  cor_wlr_dhat <- cov_wlr_dhat / (sqrt(sig2_lr) * sqrt(sig2_dhat))
  list(
    lr = lr, sig2_lr = sig2_lr, dhat = dhat_tzero,
    cov_wlr_dhat = cov_wlr_dhat, sig2_dhat = sig2_dhat, cor_wlr_dhat = cor_wlr_dhat
  )
}

#' KM difference at specified timepoints
#' @param df Data frame with survival data.
#' @param tte.name Name of time-to-event column.
#' @param event.name Name of event indicator column.
#' @param treat.name Name of treatment group column.
#' @param weight.name Name of weights column (optional).
#' @param at.points Time points for estimates.
#' @param alpha Significance level.
#' @return List with survival, difference, and CI estimates.
#' @export
KM_diff <- function(df, tte.name, event.name, treat.name, weight.name=NULL, at.points = sort(df[[tte.name]]), alpha = 0.05, seedstart = 8316951, draws = 0,
                    risk.points, draws.band = 0, tau.seq = 0.25, qtau = 0.025, show_resamples = TRUE) {

  required_cols <- c(tte.name, event.name, treat.name)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns in df:", paste(missing_cols, collapse = ", ")))
  }

  # Check treatment coding
  if (!all(df[[treat.name]] %in% c(0, 1))) {
    stop("Treatment must be numerical indicator: 0=control, 1=experimental")
  }

  # Check event coding
  if (!all(df[[event.name]] %in% c(0, 1))) {
    stop("Event must be binary (0/1).")
  }

  # Check for NA values
  if (any(is.na(df[[tte.name]]))) warning("NA values found in time-to-event column.")
  if (any(is.na(df[[event.name]]))) warning("NA values found in event column.")
  if (any(is.na(df[[treat.name]]))) warning("NA values found in treatment column.")

  tfixed <- aeqSurv(Surv(df[[tte.name]],df[[event.name]]))
  time<- tfixed[,"time"]
  delta <- tfixed[,"status"]
  z <- df[[treat.name]]
  wgt <- if (!is.null(weight.name)) df[[weight.name]] else rep(1, length(time))

  if (!all(z %in% c(0, 1))) stop("Treatment must be numerical indicator: 0=control, 1=experimental")

  if (is.unsorted(time)) {
    ord <- order(time)
    time <- time[ord]
    delta <- delta[ord]
    z <- z[ord]
    wgt <- wgt[ord]
  }

 # For simultaneous bands restrict time range
  if(draws.band > 0){
  taus <- quantile(time[delta ==1], c(qtau, 1-qtau))
  at.points<-seq(taus[1],taus[2], by =tau.seq)
  riskp <- risk.points[which(risk.points <= taus[2])]
  at.points <- sort(unique(c(at.points, riskp)))
  }

  # Treatment arm
  group_data <- extract_group_data(time, delta, wgt, z, group = 1)
  risk_event <- calculate_risk_event_counts(group_data$U, group_data$D, group_data$W, at_points = at.points, draws = draws, seedstart = seedstart)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  risk_event1 <- risk_event
  group_data1 <- group_data
  surv1 <- temp$S_KM
  sig2_surv1 <- temp$sig2_KM

  # Treatment arm
  group_data <- extract_group_data(time, delta, wgt, z, group = 0)
  risk_event <- calculate_risk_event_counts(group_data$U, group_data$D, group_data$W, at_points = at.points, draws = draws, seedstart = seedstart)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  risk_event0 <- risk_event
  group_data0 <- group_data
  surv0 <- temp$S_KM
  sig2_surv0 <- temp$sig2_KM

  dhat <- surv1 - surv0
  sig2_dhat <- sig2_surv0 + sig2_surv1

  #print(summary(sig2_dhat))

  surv0_star <- surv1_star <- dhat_star <- NULL
  c_alpha_band <- sb_lower <- sb_upper <- NULL

  if(draws.band > 0){
  # draws > 0
  set.seed(seedstart)
  # Control
  n0 <- length(group_data0$U)
  U <- group_data0$U
  W <- group_data0$W
  D <- group_data0$D
  event_mat <- outer(U, at.points, FUN = "<=")
  risk_mat  <- outer(U, at.points, FUN = ">=")
  risk_w <- colSums(risk_mat *  W)

   G0.draws <- matrix(rnorm(draws.band * n0), ncol = draws.band)
   counting_star_all <- t(event_mat * W) %*% (D * G0.draws)
   dN_star_all <- apply(counting_star_all, 2, function(x) diff(c(0, x)))
   drisk_star <- sweep(dN_star_all, 1, risk_w, "/")
   drisk_star[is.infinite(drisk_star) | is.nan(drisk_star)] <- 0
   # (length(at.points) x draws.band) dimension
   surv0_star <- (-1) * surv0 * apply(drisk_star, 2, cumsum)

  # Treatment
  n1 <- length(group_data1$U)
  U <- group_data1$U
  W <- group_data1$W
  D <- group_data1$D
  event_mat <- outer(U, at.points, FUN = "<=")
  risk_mat  <- outer(U, at.points, FUN = ">=")
  risk_w <- colSums(risk_mat *  W)
  G1.draws <- matrix(rnorm(draws.band * n1), ncol = draws.band)
  counting_star_all <- t(event_mat * W) %*% (D * G1.draws)
  dN_star_all <- apply(counting_star_all, 2, function(x) diff(c(0, x)))
  drisk_star <- sweep(dN_star_all, 1, risk_w, "/")
  drisk_star[is.infinite(drisk_star) | is.nan(drisk_star)] <- 0
  # (length(at.points) x draws.band) dimension
  surv1_star <- (-1) * surv1 * apply(drisk_star, 2, cumsum)
  dhat_star <- (surv1_star - surv0_star) / sqrt(sig2_dhat)

  #print(dim(dhat_star))

  # simultaneous band
  sups <- apply(abs(dhat_star), 2, max, na.rm = TRUE)
  c_alpha_band <- quantile(sups,c(0.95))
# Show first 20
if(show_resamples){
  matplot(at.points, dhat_star[,c(1:20)], type="s", xlab="time", ylab = "Survival differences (1st 20)", main = sprintf("c_alpha (simul. band): %.2f", c_alpha_band)
          )
}

  # simulataneous band
  sb_lower <- dhat - c_alpha_band * sqrt(sig2_dhat)
  sb_upper <- dhat + c_alpha_band * sqrt(sig2_dhat)
  }


 # Standard point-wise CIs
  c_alpha <- qnorm(1 - alpha / 2)
  lower <- dhat - c_alpha * sqrt(sig2_dhat)
  upper <- dhat + c_alpha * sqrt(sig2_dhat)

  list(
    at.points = at.points, surv0 = surv0, sig2_surv0 = sig2_surv0,
    surv1 = surv1, sig2_surv1 = sig2_surv1, dhat = dhat, sig2_dhat = sig2_dhat,
    lower = lower, upper = upper,
    dhat_star = dhat_star, surv0_star = surv0_star, surv1_star = surv1_star,
    c_alpha_band = c_alpha_band, sb_lower = sb_lower, sb_upper = sb_upper
  )
}

#' Score test statistic for survival data
#' @param nbar0 Number of events in group 0.
#' @param ybar0 Number at risk in group 0.
#' @param nbar1 Number of events in group 1.
#' @param ybar1 Number at risk in group 1.
#' @return List with z-score, score, and variance.
#' @export
z_score_calculations <- function(nbar0, ybar0, nbar1, ybar1, rho = 0, gamma = 0, S_pool = NULL){
w <- rep(1, length(nbar0))
if( rho != 0.0 | gamma != 0.0){
  if(is.null(S_pool)){
    ybar <- ybar0 + ybar1
    nbar <- nbar0 + nbar1
    dN <- diff(c(0,nbar))
    dN_Risk <- ifelse(ybar > 0, dN / ybar, 0)
    S_pool <- cumprod(1 - dN_Risk)
    S_pool <- c(1, S_pool[-length(S_pool)]) # S_pool(t-)
  }
 if(!is.null(S_pool)) S_pool <- c(1, S_pool[-length(S_pool)])
 w <- (S_pool^rho) * ((1 - S_pool)^gamma)
}

  # score test statistic
  dN.z1 <- diff(c(0, nbar1))
  dN.z0 <- diff(c(0, nbar0))
  num <- w * ybar1 * ybar0
  den <- ybar1 + ybar0
  K <- ifelse(den > 0, num / den, 0.0)
  drisk1 <- ifelse(ybar1 > 0, dN.z1 / ybar1, 0.0)
  drisk0 <- ifelse(ybar0 > 0, dN.z0 / ybar0, 0.0)
  score <- sum(K * (drisk0 - drisk1))
  i.bhat <- sum(ifelse(den > 0, (num / (den^2)) * (dN.z0 + dN.z1), 0.0))
  DNbar <- dN.z0 + dN.z1
  h1 <- ifelse(ybar1 > 0, (K^2 / ybar1), 0.0)
  h2 <- ifelse(ybar0 > 0, (K^2 / ybar0), 0.0)
  temp <- c(den - 1)
  ybar_mod <- ifelse(temp < 1, 1, temp)
  dH1 <- ifelse(ybar_mod > 0, (DNbar-1) / ybar_mod, 0.0)
  dH2 <- ifelse(den > 0, DNbar / den, 0.0)
  sig2s <- (h1+h2)*(1-dH1)*dH2
  sig2U.bzero <- sum(sig2s)
  z.score <- score / sqrt(sig2U.bzero)
  return(list(z.score=z.score, score=score, sig2.score=sig2U.bzero))
}


