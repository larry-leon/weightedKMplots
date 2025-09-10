
# ---- Utility Functions ----
#' Safe execution wrapper
#' @export
safe_run <- function(expr) {
  tryCatch(expr, error = function(e) {
    message("Error: ", e$message)
    NULL
  })
}


#' Get df_counting
#' @export
get_dfcounting <- function(df, tte.name, event.name, treat.name, arms, by.risk=12, cox.digits=3, lr.digits=3,
                           qprob=0.50, strata.name=NULL, weight.name=NULL, check.KM = TRUE, rho = 0, gamma = 0, draws = 0, seedstart = 8316951, check.seKM = TRUE) {
  safe_run({
    dfcount <- df_counting(
      df=df, tte.name=tte.name, event.name=event.name, treat.name=treat.name,
      arms=arms, by.risk=by.risk, cox.digits=cox.digits, lr.digits=lr.digits,
      qprob=qprob, strata.name=strata.name, weight.name=weight.name, check.KM = check.KM, check.seKM = check.seKM,
      rho = rho, gamma = gamma, draws = draws, seedstart=seedstart
    )
    return(dfcount)
  })
}

#' Checking results
#' @export
check_results <- function(dfcount){
  zlr_sq  <- with(dfcount,lr^2/sig2_lr)
  zCox_sq <-  with(dfcount,z.score^2)
  cat(sprintf("zlr_sq=%.6f, logrank=%.6f, zCox_sq=%.6f\n", zlr_sq, dfcount$logrank_results$chisq, zCox_sq))
}


#' Plot Kaplan-Meier curves
#' @export
plot_km <- function(df, tte.name, event.name, treat.name, weights=NULL, ...) {
  safe_run({
    surv_obj <- Surv(df[[tte.name]], df[[event.name]])
    formula <- as.formula(paste("surv_obj ~", treat.name))
    if (!is.null(weights)) {
      km_fit <- survfit(formula, data=df, weights=df[[weights]])
    } else {
      km_fit <- survfit(formula, data=df)
    }
    plot(km_fit, mark.time=TRUE, ...)
    invisible(km_fit)
  })
}

#' Plot weighted KM using custom function
#' @export
plot_weighted_km <- function(dfcount, ...) {
  safe_run({
    KM_plot_2sample_weighted_counting(
      dfcount=dfcount, risk.cex=0.725, risk_offset=0.125, risk_delta=0.05,
      show.cox=TRUE, show.logrank=FALSE, show.med=TRUE, med.font=4, ...
    )
  })
}


# Helper: Extract Group Data
#' Extract time, event, and weight data for a group
#'
#' @param time Numeric vector of times
#' @param delta Numeric vector of event indicators (1=event, 0=censored)
#' @param wgt Numeric vector of weights
#' @param z Numeric vector of group indicators
#' @param group Value of group to extract (default 1)
#' @return List with U (times), D (events), W (weights)
#' @export
extract_group_data <- function(time, delta, wgt, z, group = 1) {
  list(
    U = time[z == group],
    D = delta[z == group],
    W = wgt[z == group]
  )
}

# Helper: Calculate Risk and Event Counts
#' Calculate risk set and event counts at time points
#'
#' @param U Numeric vector of times for group
#' @param D Numeric vector of event indicators for group
#' @param W Numeric vector of weights for group
#' @param at_points Numeric vector of time points
#' @return List with ybar (risk set counts), nbar (event counts)
#' @export
calculate_risk_event_counts <- function(U, D, W, at_points, draws = 0, seedstart = 816951) {
  ybar <- colSums(outer(U, at_points, FUN = ">=") * W)
  nbar <- colSums(outer(U[D == 1], at_points, FUN = "<=") * W[D == 1])
  dN <- diff(c(0, nbar))
  # For un-weighted (all weights equal) return standard variance term
  if(length(unique(W)) == 1){
    # Greenwood
    sig2w_multiplier <- ifelse(ybar > 0 & ybar > dN, dN / (ybar * (ybar-dN)), 0.0)
    # Alternative with dN / (ybar^2)
  }
  if(length(unique(W)) > 1){
    n <- length(U)
    event_mat <- outer(U, at_points, FUN = "<=")
    risk_mat  <- outer(U, at_points, FUN = ">=")
    risk_w <- colSums(risk_mat *  W)
    result <- switch(
      as.integer(draws > 0) + 1,
      {
        # draws == 0
        counting <- colSums(event_mat * (D * W))
        dN_w <- diff(c(0, counting))
        dJ <- ifelse(risk_w == 1, 0, (dN_w - 1) / (risk_w - 1))
        dL <- ifelse(risk_w == 0, 0, dN_w / risk_w)
        h2 <- ifelse(risk_w > 0, (1 / (risk_w)), 0)
        sig2w_multiplier <- (h2 * (dN_w - dL))^2
        sig2w_multiplier
      },
      {
        # draws > 0
        set.seed(seedstart)
        G.draws <- matrix(rnorm(draws * n), ncol = draws)
        counting_star_all <- t(event_mat * W) %*% (D * G.draws)
        dN_star_all <- apply(counting_star_all, 2, function(x) diff(c(0, x)))
        drisk_star <- sweep(dN_star_all, 1, risk_w, "/")
        drisk_star[is.infinite(drisk_star) | is.nan(drisk_star)] <- 0
        sig2w_multiplier <- apply(drisk_star, 1, var)
        sig2w_multiplier
      }
    )
   }
  list(ybar = ybar, nbar = nbar, sig2w_multiplier = sig2w_multiplier)
}

# Helper : Get Censoring and Event Times
#' Get censoring and event times and their indices
#'
#' @param time Numeric vector of times
#' @param delta Numeric vector of event indicators
#' @param z Numeric vector of group indicators
#' @param group Value of group to extract
#' @param censoring_allmarks Logical; if FALSE, remove events from censored
#' @param at_points Numeric vector of time points
#' @return List with cens (censored times), ev (event times), idx_cens, idx_ev, idx_ev_full
#' @export
get_censoring_and_events <- function(time, delta, z, group, censoring_allmarks, at_points) {
  cens <- time[z == group & delta == 0]
  ev <- sort(unique(time[z == group & delta == 1]))
  if (!censoring_allmarks) cens <- setdiff(cens, ev)
  idx_cens <- match(cens, at_points)
  idx_ev <- match(ev, at_points)
  ev <- c(ev, max(time[z == group]))
  idx_ev_full <- match(ev, at_points)
  list(
    cens = cens,
    ev = ev,
    idx_cens = idx_cens,
    idx_ev = idx_ev,
    idx_ev_full = idx_ev_full
  )
}

# Helper : Get Risk Points
#' Get risk set counts at specified risk points
#'
#' @param ybar Numeric vector of risk set counts
#' @param risk_points Numeric vector of risk points
#' @param at_points Numeric vector of time points
#' @return Numeric vector of risk set counts at risk points
#' @export
get_riskpoints <- function(ybar, risk_points, at_points) {
  ybar[match(risk_points, at_points)]
}


#' Creates a counting process dataset for survival analysis
#'
#' This function prepares a dataset for survival analysis using the counting process approach,
#' including risk set, event counts, Kaplan-Meier estimates, log-rank and Cox model results,
#' and quantile estimates for two groups (treatment and control).
#'
#' @param df Data frame containing the survival data.
#' @param tte.name Name of the time-to-event variable (string).
#' @param event.name Name of the event indicator variable (string, 1=event, 0=censored).
#' @param treat.name Name of the treatment group variable (string, 0=control, 1=treatment).
#' @param weight.name Optional name of the weights variable (string, default NULL).
#' @param strata.name Optional name of the stratification variable (string, default NULL).
#' @param arms Character vector of length 2 with names for treatment and control arms.
#' @param time.zero Time value to use as zero (default 0).
#' @param tpoints.add Additional time points to include (numeric vector, default 0).
#' @param by.risk Interval for risk table (default 6).
#' @param time.zero.label Label for time zero (default 0.0).
#' @param risk.add Additional risk points (numeric vector, default NULL).
#' @param get.cox Logical; whether to compute Cox model results (default TRUE).
#' @param cox.digits Number of digits for Cox model results (default 2).
#' @param lr.digits Number of digits for log-rank results (default 2).
#' @param cox.eps Threshold for Cox p-value formatting (default 0.001).
#' @param lr.eps Threshold for log-rank p-value formatting (default 0.001).
#' @param qprob Quantile probability for median/quantile estimation (default 0.5).
#' @param rho Weighting parameter for log-rank test (default 0).
#' @param gamma Weighting parameter for log-rank test (default 0).
#' @param conf_level Confidence level for quantile CI (default 0.95).
#' @param check.KM Logical; whether to check KM curve fits (default TRUE).
#' @param stop.onerror Logical; whether to stop on error (default FALSE).
#' @param censoring_allmarks Logical; whether to mark all censoring times (default TRUE).
#' @return A list containing risk set, event counts, KM estimates, log-rank and Cox results, quantiles, and more.
#' @export
df_counting <- function(df, tte.name, event.name, treat.name, weight.name=NULL, strata.name = NULL, arms=c("treat","control"),
                        time.zero=0, tpoints.add=c(0),
                        by.risk=6, time.zero.label = 0.0, risk.add=NULL, get.cox=TRUE, cox.digits=2, lr.digits=2,
                        cox.eps = 0.001, lr.eps = 0.001, verbose = FALSE,
                        qprob=0.5, rho = 0, gamma = 0, scheme = "fh",
                        conf_level = 0.95, check.KM = TRUE, check.seKM = TRUE, draws = 0, seedstart = 8316951,
                        stop.onerror=FALSE,censoring_allmarks=TRUE) {

  validate_input(df, c(tte.name, event.name, treat.name, weight.name))

  ans <- list()

  cox_results <- NULL
  if (get.cox) {
    if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required for Cox model.")
    treat.name.strata <- treat.name
    if (!is.null(strata.name)) {
      treat.name.strata <- paste(treat.name,"+",paste("strata(",eval(strata.name),")"))
    }
    cox_formula <- as.formula(paste0("survival::Surv(", tte.name, ",", event.name, ") ~ ", treat.name.strata))
    if (!is.null(weight.name)) {
      cox_fit <- survival::coxph(cox_formula, data = df, weights = df[[weight.name]], robust = TRUE)
    } else {
      cox_fit <- survival::coxph(cox_formula, data = df, robust=TRUE)
    }
    cox_summary <- summary(cox_fit)
    cox_score <- cox_summary$sctest[1]
    hr <- exp(cox_fit$coef)
    hr_ci <- exp(confint(cox_fit))
    pval <- cox_summary$coefficients[1, "Pr(>|z|)"]
    cox_text <- paste0("HR = ", round(hr, cox.digits),
                       " (", round(hr_ci[1], cox.digits), ", ", round(hr_ci[2], cox.digits), ")",
                       ", p = ", format_pval(pval, eps = cox.eps, digits = cox.digits))
    cox_results <- list(
      cox_fit = cox_fit,
      hr = hr,
      hr_ci = hr_ci,
      pval = pval,
      score = cox_score,
      cox_text = cox_text
    )
  }
  ans$cox_results <- cox_results
  logrank_results <- NULL
  if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required for logrank test.")
  surv_obj <- survival::Surv(df[[tte.name]], df[[event.name]])
  group <- df[[treat.name]]
 # survdiff does not support case weights
  if (!is.null(strata.name)) {
    strata_var <- df[[strata.name]]
    logrank_fit <- survival::survdiff(surv_obj ~ group + strata(strata_var), rho = rho)
  } else {
    logrank_fit <- survival::survdiff(surv_obj ~ group, rho = rho)
  }

  chisq <- logrank_fit$chisq
  pval <- 1 - pchisq(chisq, df = 1)
  logrank_text <- paste0("Logrank p = ", format_pval(pval, eps = lr.eps, digits = lr.digits))
  logrank_results <- list(
    chisq = chisq,
    pval = pval,
    logrank_text = logrank_text
  )

  ans$logrank_results <- logrank_results

  # Implement timefix per survival package
  tfixed <- aeqSurv(Surv(df[[tte.name]],df[[event.name]]))
  time<- tfixed[,"time"]
  delta <- tfixed[,"status"]
  z <- df[[treat.name]]

  strata <- if (!is.null(strata.name)) df[[strata.name]] else rep("All", length(time))
  wgt <- if (!is.null(weight.name)) df[[weight.name]] else rep(1, length(time))

 if (!all(z %in% c(0, 1))) stop("Treatment must be numerical indicator: 0=control, 1=experimental")

  if (is.unsorted(time)) {
    ord <- order(time)
    time <- time[ord]
    delta <- delta[ord]
    z <- z[ord]
    strata <- strata[ord]
    wgt <- wgt[ord]
  }
  stratum <- unique(strata)
  risk.points <- sort(unique(c(seq(time.zero.label, max(time), by = ifelse(is.null(by.risk), 1, by.risk)), risk.add)))
  ans$risk.points <- risk.points
  ans$risk.points.label <- as.character(c(time.zero.label, risk.points[-1]))

  at_points <- sort(unique(c(time, time.zero, tpoints.add, risk.points)))

  # Treatment arm
  group_data <- extract_group_data(time, delta, wgt, z, group = 1)
  risk_event <- calculate_risk_event_counts(group_data$U, group_data$D, group_data$W, at_points, draws, seedstart)
  cens_ev <- get_censoring_and_events(time, delta, z, 1, censoring_allmarks, at_points)
  riskpoints1 <- get_riskpoints(risk_event$ybar, risk.points, at_points)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  surv1 <- temp$S_KM
  sig2_surv1 <- temp$sig2_KM

  nbar1 <- risk_event$nbar
  ybar1 <- risk_event$ybar

  # Store in ans
  ans$idx1 <- cens_ev$idx_cens
  ans$idv1 <- cens_ev$idx_ev_full
  idv1.check <- cens_ev$idx_ev
  ans$idv1.check <- idv1.check
  ans$ev1 <- cens_ev$ev
  ans$cens1 <- cens_ev$cens
  ans$riskpoints1 <- riskpoints1

  # Control arm
  group_data <- extract_group_data(time, delta, wgt, z, group = 0)
  risk_event <- calculate_risk_event_counts(group_data$U, group_data$D, group_data$W, at_points, draws, seedstart)
  cens_ev <- get_censoring_and_events(time, delta, z, 0, censoring_allmarks, at_points)
  riskpoints0 <- get_riskpoints(risk_event$ybar, risk.points, at_points)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  surv0 <- temp$S_KM
  sig2_surv0 <- temp$sig2_KM

  nbar0 <- risk_event$nbar
  ybar0 <- risk_event$ybar

  # Store in ans
  ans$idx0 <- cens_ev$idx_cens
  ans$idv0 <- cens_ev$idx_ev_full
  idv0.check <- cens_ev$idx_ev
  ans$idv0.check <- idv0.check
  ans$ev0 <- cens_ev$ev
  ans$cens0 <- cens_ev$cens
  ans$riskpoints0 <- riskpoints0

  # Pooled KM estimates
  risk_event <- calculate_risk_event_counts(time, delta, wgt, at_points)
  temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
  survP <- temp$S_KM
  sig2_survP <- temp$sig2_KM
  rm("temp")

  get_lr <- wlr_estimates(ybar0 = ybar0, ybar1 = ybar1, nbar0 = nbar0, nbar1 = nbar1,
                          rho = rho, gamma = gamma, S_pool = survP, scheme = scheme)
  # Quantiles
  get_kmq <- km_quantile_table(at_points, surv0, se0 = sqrt(sig2_surv0), surv1, se1 = sqrt(sig2_surv1), arms,
                               qprob = qprob, type = c("midpoint"), conf_level = conf_level)
  ans$quantile_results <- get_kmq

  get_score <- z_score_calculations(nbar0, ybar0, nbar1, ybar1, rho = rho, gamma = gamma, S_pool = survP, scheme = scheme)
  ans$z.score <- get_score$z.score

  # KM curve checks
  if (check.KM) {
    # control arm will be first stratum
    check_fit <- survfit(Surv(time,delta) ~ z, weights = wgt)
    check_sfit <- summary(check_fit)
    strata_names <- as.character(check_sfit$strata)
    strata_lengths <- rle(strata_names)$lengths
    strata_labels <- rle(strata_names)$values
    split_times <- split(check_sfit$time, rep(strata_labels, strata_lengths))
    split_surv  <- split(check_sfit$surv, rep(strata_labels, strata_lengths))
    split_se  <- split(check_sfit$std.err, rep(strata_labels, strata_lengths))
    df0_check <- data.frame(time=split_times[[1]], surv=split_surv[[1]], se = split_se[[1]])
    df1_check <- data.frame(time=split_times[[2]], surv=split_surv[[2]], se = split_se[[2]])
    # Note: the quantile() function can yield different results than median table
    # The median table calculations appear more stable ...
    if(qprob != 0.50){
      qtab <- quantile(check_fit,probs=c(qprob))
      strata_names <- rownames(qtab$quantile)
      quantile_table <- data.frame(
        stratum = rep(strata_names, each = length(qprob)),
        quantile = rep(qprob, times = length(strata_names)),
        time = as.vector(qtab$quantile),
        lower = as.vector(qtab$lower),
        upper = as.vector(qtab$upper)
      )
    ans$quantile_check <- quantile_table
    }
    if(qprob == 0.50){
      qtab <- summary(check_fit)$table[, c("median", "0.95LCL", "0.95UCL")]
      quantile_table <- data.frame(
        time = qtab[, "median"],
        lower = qtab[, "0.95LCL"],
        upper = qtab[, "0.95UCL"]
      )
      ans$quantile_check <- quantile_table
        }
    # First row is control here
    qcheck_0 <- quantile_table[1,c("time","lower","upper")]
    aa <- c(unlist(qcheck_0))
    bb <- c(unlist(get_kmq[2,c("quantile","lower","upper")]))
    dcheck <- round(abs(aa-bb),6)
    if(max(dcheck,na.rm=TRUE) > 1e-6){
      msg <- paste0(arms[2]," : ", "Control: discrepancy in quantile calculations")
    if(verbose){  if (stop.onerror) stop(msg) else warning(msg)
    }
    }
    qcheck_1 <- quantile_table[2,c("time","lower","upper")]
    aa <- c(unlist(qcheck_1))
    bb <- c(unlist(get_kmq[1,c("quantile","lower","upper")]))
    dcheck <- round(abs(aa-bb),6)
    if(max(dcheck,na.rm=TRUE) > 1e-6){
      msg <- paste0(arms[1]," : ", "Treatment: discrepancy in quantile calculations")
    if(verbose){  if (stop.onerror) stop(msg) else warning(msg)
    }
    }
    check_km_curve <- function(time, S.KM, se.KM, df_check, group_name = "Group", check.seKM = TRUE) {
        if (any(S.KM < 0 | S.KM > 1)) {
        msg <- paste0(group_name, " : ","KM curve has values outside [0,1].")
        if (stop.onerror) stop(msg) else warning(msg)
      }
      if (any(diff(S.KM) > 0)) {
        msg <- paste0(group_name, " : ", " KM curve is not non-increasing.")
        if (stop.onerror) stop(msg) else warning(msg)
      }

      if(round(max(abs(S.KM-df_check$surv)),8)) {
        msg <- paste0(group_name," : ", "Discrepancy in KM curve fit.")
        if (stop.onerror) stop(msg) else warning(msg)
      }

      if(round(max(abs(se.KM-df_check$se)),8) && check.seKM) {
        yymax <- max(c(se.KM, df_check$se))
        plot(time,se.KM, type="s", lty=1, col="lightgrey", lwd=4, ylim=c(0,yymax), xlab="time", ylab="SE(KM)")
        with(df_check, lines(time, se, type="s", lty=2, lwd=1, col="red"))
        legend("topleft",c("Mine","Survfit"), lty=c(1,2),col=c("lightgrey","red"), lwd=c(4,1),bty="n", cex=0.8)
        title(main=group_name)


      if(verbose){  msg <- paste0(group_name," : ", "Discrepancy in se(KM) curve fit.")
      ratio <- se.KM / with(df_check,se)
      print(summary(ratio))
      if (stop.onerror) stop(msg) else warning(msg)
      }
      }
    }
    par(mfrow=c(1,2))
    check_km_curve(at_points[idv0.check],surv0[idv0.check], sqrt(sig2_surv0[idv0.check]), df0_check, "control", check.seKM = check.seKM)
    check_km_curve(at_points[idv1.check],surv1[idv1.check], sqrt(sig2_surv1[idv1.check]), df1_check, "treat", check.seKM = check.seKM)
  }
  ans$lr <- get_lr$lr
  ans$sig2_lr <- get_lr$sig2
  ans$at.points <- at_points
  ans$strata <- strata
  ans$ybar0 <- ybar0
  ans$nbar0 <- nbar0
  ans$surv0 <- surv0
  ans$sig2_surv0 <- sig2_surv0
  ans$ybar1 <- ybar1
  ans$nbar1 <- nbar1
  ans$surv1 <- surv1
  ans$sig2_surv1 <- sig2_surv1
  ans$survP <- survP
  ans$sig2_survP <- sig2_survP

  if (length(stratum) > 1) {
    ybar0_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    nbar0_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    ybar1_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    nbar1_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    surv0_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    surv1_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    sig2_surv0_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    sig2_surv1_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    survP_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    sig2_survP_mat <- matrix(NA, nrow = length(at_points), ncol = length(stratum))
    lr_stratified <- 0.0
    sig2_lr_stratified <- 0.0
    score_stratified <- 0.0
    sig2_score_stratified <- 0.0
    for (ss in seq_along(stratum)) {
      this_stratum <- stratum[ss]
      U0_s <- time[z == 0 & strata == this_stratum]
      D0_s <- delta[z == 0 & strata == this_stratum]
      W0_s <- wgt[z == 0 & strata == this_stratum]

      U1_s <- time[z == 1 & strata == this_stratum]
      D1_s <- delta[z == 1 & strata == this_stratum]
      W1_s <- wgt[z == 1 & strata == this_stratum]

      risk_event <- calculate_risk_event_counts(U1_s, D1_s, W1_s, at_points, draws, seedstart)
      temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
      surv1_mat[, ss] <- temp$S_KM
      sig2_surv1_mat[, ss] <- temp$sig2_KM
      rm("temp")

      nbar1_s <- risk_event$nbar
      ybar1_s <- risk_event$ybar

      risk_event <- calculate_risk_event_counts(U0_s, D0_s, W0_s, at_points, draws, seedstart)
      temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
      surv0_mat[, ss] <- temp$S_KM
      sig2_surv0_mat[, ss] <- temp$sig2_KM
      rm("temp")

      nbar0_s <- risk_event$nbar
      ybar0_s <- risk_event$ybar

      U_s <- time[strata == this_stratum]
      D_s <- delta[strata == this_stratum]
      W_s <- wgt[strata == this_stratum]

      risk_event <- calculate_risk_event_counts(U_s, D_s, W_s, at_points, draws, seedstart)
      temp <- KM_estimates(ybar = risk_event$ybar, nbar = risk_event$nbar, sig2w_multiplier = risk_event$sig2w_multiplier)
      survP_mat[, ss] <- temp$S_KM
      sig2_survP_mat[, ss] <- temp$sig2_KM
      S_pool <- temp$S_KM

      get_score <- z_score_calculations(nbar0 = nbar0_s,ybar0 = ybar0_s,nbar1 = nbar1_s, ybar1 = ybar1_s,
                                        rho = rho, gamma = gamma, S_pool = S_pool)

      score_stratified <- score_stratified + get_score$score
      sig2_score_stratified <- sig2_score_stratified + get_score$sig2.score

      temp <- wlr_estimates(ybar0 = ybar0_s, ybar1 = ybar1_s, nbar0 = nbar0_s, nbar1 = nbar1_s,
                            rho = rho, gamma = gamma, S_pool = S_pool)
      lr_stratified <- lr_stratified + temp$lr
      sig2_lr_stratified <- sig2_lr_stratified + temp$sig2

      ybar0_mat[, ss] <- ybar0_s
      nbar0_mat[, ss] <- nbar0_s
      ybar1_mat[, ss] <- ybar1_s
      nbar1_mat[, ss] <- nbar1_s
 }

    # For each observation, get the index of its time in at_points
    id_time <- match(time, at_points)  # 'time' is the vector of observed times
    # For each observation, get the index of its stratum in stratum
    id_stratum <- vapply(strata, function(x) which(stratum == x), integer(1))
    # Now extract for each observation
    ans$ybar0_stratum      <- ybar0_mat[cbind(id_time, id_stratum)]
    ans$nbar0_stratum      <- nbar0_mat[cbind(id_time, id_stratum)]
    ans$ybar1_stratum      <- ybar1_mat[cbind(id_time, id_stratum)]
    ans$nbar1_stratum      <- nbar1_mat[cbind(id_time, id_stratum)]
    ans$surv0_stratum      <- surv0_mat[cbind(id_time, id_stratum)]
    ans$sig2_surv0_stratum <- sig2_surv0_mat[cbind(id_time, id_stratum)]
    ans$surv1_stratum      <- surv1_mat[cbind(id_time, id_stratum)]
    ans$sig2_surv1_stratum <- sig2_surv1_mat[cbind(id_time, id_stratum)]
    ans$survP_stratum      <- survP_mat[cbind(id_time, id_stratum)]
    ans$sig2_survP_stratum <- sig2_survP_mat[cbind(id_time, id_stratum)]

    ans$lr_stratified <- lr_stratified
    ans$sig2_lr_stratified <- sig2_lr_stratified
    ans$z.score_stratified <- score_stratified/sqrt(sig2_score_stratified)

    }

  # Compare with survdiff for gamma=0 (survdiff only handles gamma=0)
  if(is.null(weight.name) && is.null(strata.name) && gamma == 0){
    z_lr <- with(get_lr,lr/sqrt(sig2))
    zsq_lr_check <- logrank_results$chisq
    if(round(z_lr^2 - zsq_lr_check,8)>0){
      warning("Discrepancy with log-rank and survdiff")
      cat("Discrepancy with survdiff",c(z_lr^2,zsq_lr_check),"\n")
    }
  }

  if(is.null(weight.name) && !is.null(strata.name) && gamma == 0){
    z_lr <- lr_stratified/sqrt(sig2_lr_stratified)
    zsq_lr_check <- logrank_results$chisq
    if(round(z_lr^2 - zsq_lr_check,8)>0){
      warning("Discrepancy with log-rank and survdiff")
      cat("Discrepancy with survdiff",c(z_lr^2,zsq_lr_check),"\n")
    }
  }

  ans
}


