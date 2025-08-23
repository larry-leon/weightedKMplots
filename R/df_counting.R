
#' Creates counting process dataset
#' @export
df_counting <- function(df, tte.name, event.name, treat.name, weight.name=NULL, strata.name = NULL, arms=c("treat","control"), time.zero=0, tpoints.add=c(0),
                        by.risk=6, time.zero.label = 0.0, risk.add=NULL, get.cox=TRUE, cox.digits=2, lr.digits=2, cox.eps = 0.001, lr.eps = 0.001,
                        qprob=0.5, conf_level = 0.95, check.KM=TRUE,stop.onerror=FALSE,censoring_allmarks=TRUE) {

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
      cox_fit <- survival::coxph(cox_formula, data = df, weights = df[[weight.name]],robust=TRUE)
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
  if (!is.null(strata.name)) {
    strata_var <- df[[strata.name]]
    logrank_formula <- as.formula(paste0("surv_obj ~ group + strata(strata_var)"))
    logrank_fit <- eval(bquote(survival::survdiff(.(logrank_formula))))
  } else {
    logrank_fit <- survival::survdiff(surv_obj ~ group)
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

  U0 <- time[z == 0]
  D0 <- delta[z == 0]
  W0 <- wgt[z == 0]

  ybar0 <- colSums(outer(U0, at_points, FUN = ">=") * W0)
  nbar0 <- colSums(outer(U0[D0 == 1], at_points, FUN = "<=") * W0[D0 == 1])

  cens0 <- time[z==0 & delta == 0]
  ev0 <- sort(unique(time[z==0 & delta == 1]))
  # Censoring that are NOT an event
  if(!censoring_allmarks) cens0 <- setdiff(cens0, ev0)
  idx0 <- match(cens0, at_points)
  ans$idx0 <- idx0
  # used for checking KM fits
  idv0.check <- match(ev0, at_points)
  # append events to include max time
  ev0 <- c(ev0,max(time[z==0]))
  idv0 <- match(ev0, at_points)
  ans$idv0 <- idv0
  ans$ev0 <- ev0
  ans$cens0 <- cens0
  ans$riskpoints0 <- ybar0[match(risk.points,at_points)]


  temp <- KM_estimates(ybar = ybar0, nbar = nbar0)
  surv0 <- temp$S_KM
  sig2_surv0 <- temp$sig2_KM
  U1 <- time[z == 1]
  D1 <- delta[z == 1]
  W1 <- wgt[z== 1]

  ybar1 <- colSums(outer(U1, at_points, FUN = ">=") * W1)
  nbar1 <- colSums(outer(U1[D1 == 1], at_points, FUN = "<=") * W1[D1 == 1])

  cens1 <- time[z==1 & delta == 0]
  ev1 <- sort(unique(time[z==1 & delta == 1]))
  if(!censoring_allmarks) cens1 <- setdiff(cens1, ev1)
  idx1 <- match(cens1, at_points)
  ans$idx1 <- idx1
  idv1.check <- match(ev1, at_points)
  ev1 <- c(ev1,max(time[z==1]))
  idv1 <- match(ev1, at_points)
  ans$idv1 <- idv1

  ans$ev1 <- ev1
  ans$cens1 <- cens1
  ans$riskpoints1 <- ybar1[match(risk.points,at_points)]

  temp <- KM_estimates(ybar = ybar1, nbar = nbar1)
  surv1 <- temp$S_KM
  sig2_surv1 <- temp$sig2_KM

  # Pooled KM estimates
  temp <- KM_estimates(ybar = ybar0 + ybar1, nbar = nbar0 + nbar1)
  survP <- temp$S_KM
  sig2_survP <- temp$sig2_KM
  get_lr <- wlr_estimates(ybar0 = ybar0, ybar1 = ybar1, nbar0 = nbar0, nbar1 = nbar1, rho = 0, gamma = 0)
  # Quantiles
  get_kmq <- km_quantile_table(at_points, surv0, se0=sqrt(sig2_surv0), surv1, se1=sqrt(sig2_surv1), arms, qprob = qprob, type = c("midpoint"), conf_level = conf_level)
  ans$quantile_results <- get_kmq

  get_score <- z_score_calculations(nbar0,ybar0,nbar1,ybar1)
  ans$z.score <- get_score$z.score

  # KM curve checks
  if (check.KM) {
    # control arm will be first stratum
    check_fit <- survfit(Surv(time,delta) ~ z, weights=wgt)
    check_sfit <- summary(check_fit)
    strata_names <- as.character(check_sfit$strata)
    strata_lengths <- rle(strata_names)$lengths
    strata_labels <- rle(strata_names)$values
    split_times <- split(check_sfit$time, rep(strata_labels, strata_lengths))
    split_surv  <- split(check_sfit$surv, rep(strata_labels, strata_lengths))
    df0_check <- data.frame(time=split_times[[1]], surv=split_surv[[1]])
    df1_check <- data.frame(time=split_times[[2]], surv=split_surv[[2]])
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
    }
    if(qprob == 0.50){
      qtab <- summary(check_fit)$table[, c("median", "0.95LCL", "0.95UCL")]
      quantile_table <- data.frame(
        time = qtab[, "median"],
        lower = qtab[, "0.95LCL"],
        upper = qtab[, "0.95UCL"]
      )
    }
    # First row is control here
    qcheck_0 <- quantile_table[1,c("time","lower","upper")]
    aa <- c(unlist(qcheck_0))
    bb <- c(unlist(get_kmq[2,c("quantile","lower","upper")]))
    dcheck <- round(abs(aa-bb),6)
    if(max(dcheck,na.rm=TRUE) > 1e-6){
      msg <- paste0(arms[2]," : ", "discrepancy in quantile calculations")
      if (stop.onerror) stop(msg) else warning(msg)
    }
    qcheck_1 <- quantile_table[2,c("time","lower","upper")]
    aa <- c(unlist(qcheck_1))
    bb <- c(unlist(get_kmq[1,c("quantile","lower","upper")]))
    dcheck <- round(abs(aa-bb),6)
    if(max(dcheck,na.rm=TRUE) > 1e-6){
      msg <- paste0(arms[1]," : ", "discrepancy in quantile calculations")
      if (stop.onerror) stop(msg) else warning(msg)
    }
    check_km_curve <- function(S.KM, df_check, group_name = "Group") {
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
    }
    check_km_curve(surv1[idv1.check], df1_check, "Group 1")
    check_km_curve(surv0[idv0.check], df0_check, "Group 0")
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

      ybar0_s <- colSums(outer(U0_s, at_points, FUN = ">=") * W0_s)
      nbar0_s <- colSums(outer(U0_s[D0_s == 1], at_points, FUN = "<=") * W0_s[D0_s == 1])

      ybar1_s <- colSums(outer(U1_s, at_points, FUN = ">=") * W1_s)
      nbar1_s <- colSums(outer(U1_s[D1_s == 1], at_points, FUN = "<=") * W1_s[D1_s == 1])

      get_score <- z_score_calculations(nbar0 = nbar0_s,ybar0 = ybar0_s,nbar1 = nbar1_s, ybar1 = ybar1_s)
      score_stratified <- score_stratified + get_score$score
      sig2_score_stratified <- sig2_score_stratified + get_score$sig2.score

      ybar0_mat[, ss] <- ybar0_s
      nbar0_mat[, ss] <- nbar0_s
      ybar1_mat[, ss] <- ybar1_s
      nbar1_mat[, ss] <- nbar1_s
      temp <- KM_estimates(ybar = ybar0_s, nbar = nbar0_s)
      surv0_mat[, ss] <- temp$S_KM
      sig2_surv0_mat[, ss] <- temp$sig2_KM
      temp <- KM_estimates(ybar = ybar1_s, nbar = nbar1_s)
      surv1_mat[, ss] <- temp$S_KM
      sig2_surv1_mat[, ss] <- temp$sig2_KM
      temp <- KM_estimates(ybar = ybar0_s + ybar1_s, nbar = nbar0_s + nbar1_s)
      survP_mat[, ss] <- temp$S_KM
      sig2_survP_mat[, ss] <- temp$sig2_KM
      temp <- wlr_estimates(ybar0 = ybar0_s, ybar1 = ybar1_s, nbar0 = nbar0_s, nbar1 = nbar1_s, rho = 0, gamma = 0)
      lr_stratified <- lr_stratified + temp$lr
      sig2_lr_stratified <- sig2_lr_stratified + temp$sig2
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

  if(is.null(weight.name) && is.null(strata.name)){
    z_lr <- with(get_lr,lr/sqrt(sig2))
    zsq_lr_check <- logrank_results$chisq
    if(round(z_lr^2 - zsq_lr_check,8)>0){
      warning("Discrepancy with log-rank and survdiff")
      cat("Discrepancy with survdiff",c(z_lr^2,zsq_lr_check),"\n")
    }
  }

  if(is.null(weight.name) && !is.null(strata.name)){
    z_lr <- lr_stratified/sqrt(sig2_lr_stratified)
    zsq_lr_check <- logrank_results$chisq
    if(round(z_lr^2 - zsq_lr_check,8)>0){
      warning("Discrepancy with log-rank and survdiff")
      cat("Discrepancy with survdiff",c(z_lr^2,zsq_lr_check),"\n")
    }
  }

  ans
}


