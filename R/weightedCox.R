
get_cphw <- function(Y,Delta,X,wt,alpha=0.025){
  calpha <- qnorm(1-alpha)
  toget <- which(wt>0)
  yy <- Y[toget]
  dd <- Delta[toget]
  xx <- X[toget]
  ww <- wt[toget]
  fit  <- coxph(Surv(yy,dd)~xx + offset(-log(ww)), weights=c(ww))
  # Weights corresponding to events
  #w <- coxph.detail(fit)$nevent.wt
  w <- coxph.detail(fit)$nevent
  bhat <- fit$coefficients
   # For variance
  temp <- coxph(Surv(yy,dd)~ xx, init=c(fit$coefficients), iter=0)
  temp2 <- coxph.detail(temp)
  v <- temp2$imat
  A <- (sum(w*v))^-1
  B <- sum(v*(w^2))
  sebhat <- sqrt(A*B*A)
  lower <- bhat - calpha*sebhat
  upper <- bhat + calpha*sebhat
  return(list(bhat=bhat,sebhat=sebhat,lower=lower,upper=upper))
}


# Redefine cph.weighted with MB/custom_time fix
cph.weighted <- function(
    df,
    scheme = "fh",
    scheme_params = list(rho = 0, gamma = 0.5),
    outcome.name = "tte",
    event.name = "event",
    treat.name = "treat",
    timefix = TRUE,
    return_model = FALSE
) {
  supported_schemes <- c("fh", "schemper", "XO", "MB", "custom_time", "fh_exp1", "fh_exp2")
  if (!(scheme %in% supported_schemes)) {
    stop("scheme must be one of: ", paste(supported_schemes, collapse = ", "))
  }
  required_cols <- c(outcome.name, event.name, treat.name)
  if (!all(required_cols %in% colnames(df))) {
    stop("Data frame must contain columns: ", paste(required_cols, collapse = ", "))
  }
  Y <- df[[outcome.name]]
  Delta <- df[[event.name]]
  X <- df[[treat.name]]
  id0 <- order(Y)
  Y <- Y[id0]; Delta <- Delta[id0]; X <- X[id0]
  if (timefix) {
    tfixed <- aeqSurv(Surv(Y, Delta))
    Y <- tfixed[,"time"]
    Delta <- tfixed[,"status"]
  }
  atpoints <- Y
  km.pool <- survival::survfit(survival::Surv(Y, Delta) ~ 1)
  S.pool <- summary(km.pool, times = atpoints)$surv
  # Only include tpoints in main list for schemes that do NOT require it in scheme_params
  if (scheme %in% c("MB", "custom_time")) {
    scheme_params$tpoints <- atpoints
    wt_args <- c(list(S = S.pool, scheme = scheme), scheme_params)
  } else {
    wt_args <- c(list(S = S.pool, scheme = scheme, tpoints = atpoints), scheme_params)
  }
  # Parameter validation for each scheme
  if (scheme == "fh" && (is.null(scheme_params$rho) || is.null(scheme_params$gamma))) {
    stop("For Fleming-Harrington weights, specify both rho and gamma in scheme_params.")
  }
  if (scheme == "schemper" && (is.null(scheme_params$Scensor) || length(scheme_params$Scensor) != length(S.pool))) {
    stop("For Schemper weights, provide Scensor (censoring KM) of same length as S in scheme_params.")
  }
  if (scheme == "XO" && (is.null(scheme_params$Ybar) || length(scheme_params$Ybar) != length(S.pool))) {
    stop("For XO weights, provide Ybar (risk set sizes) of same length as S in scheme_params.")
  }
  if (scheme == "MB" && (is.null(scheme_params$tpoints) || length(scheme_params$tpoints) != length(S.pool) || is.null(scheme_params$mb_tstar))) {
    stop("For MB weights, provide tpoints (time points) of same length as S and mb_tstar (cutoff time) in scheme_params.")
  }
  if (scheme == "custom_time" && (is.null(scheme_params$tpoints) || length(scheme_params$tpoints) != length(S.pool) || is.null(scheme_params$t.tau))) {
    stop("For custom_time weights, provide tpoints (time points) of same length as S and t.tau (cutoff time) in scheme_params.")
  }
  # Calculate weights
  wt <- tryCatch(
    do.call(wt.rg.S, wt_args),
    error = function(e) {
      stop(sprintf("Weight calculation failed for scheme %s: %s", scheme, e$message))
    }
  )
  # Fit weighted Cox model
  fit_wt <- tryCatch(
    get_cphw(Y, Delta, X, wt),
    error = function(e) {
      stop(sprintf("Cox model failed for scheme %s: %s", scheme, e$message))
    }
  )
  result <- data.frame(
    scheme = scheme,
    hr = exp(fit_wt$bhat),
    var = (exp(fit_wt$bhat)^2) * (fit_wt$sebhat^2),
    lower = exp(fit_wt$lower),
    upper = exp(fit_wt$upper)
  )
  if (return_model) attr(result, "model") <- fit_wt
  return(result)
}




cox_score_rhogamma <- function(beta, time, delta, z, w_hat = rep(1,length(time)), wt_rg = rep(1,length(time))) {
  at_points <- time
  tt0 <- time[z == 0]
  dd0 <- delta[z == 0]
  w0_hat <- w_hat[z == 0]
  risk_z0 <- colSums(outer(tt0, at_points, FUN = ">=") * w0_hat)
  event_mat0 <- outer(tt0[dd0 == 1], at_points, FUN = "<=") * w0_hat[dd0 == 1]
  counting0 <- colSums(event_mat0)
  dN_z0 <- diff(c(0, counting0))
  tt1 <- time[z == 1]
  dd1 <- delta[z == 1]
  w1_hat <- w_hat[z == 1]
  risk_z1 <- colSums(outer(tt1, at_points, FUN = ">=") * w1_hat * exp(beta))
  event_mat1 <- outer(tt1[dd1 == 1], at_points, FUN = "<=") * w1_hat[dd1 == 1]
  counting1 <- colSums(event_mat1)
  dN_z1 <- diff(c(0, counting1))
  num <- wt_rg * risk_z1 * risk_z0
  den <- risk_z0 + risk_z1
  K <- ifelse(den > 0, num / den, 0.0)
  drisk1 <- ifelse(risk_z1 > 0, dN_z1 / risk_z1, 0.0)
  drisk0 <- ifelse(risk_z0 > 0, dN_z0 / risk_z0, 0.0)
  score <- sum(K * (drisk0 - drisk1))
  return(score)
}


cox_rhogamma_old <- function(dfcount, scheme = "fh", scheme_params = list(rho = 0, gamma = 0.5)) {
  time <- dfcount$time
  delta <- dfcount$delta
  z <- dfcount$z
  w_hat <- dfcount$w_hat
  atpoints <- time

  S.pool <- dfcount$survP_all
  G.pool <- dfcount$survG_all

  stopifnot(is.numeric(time), is.numeric(delta), is.numeric(z), is.numeric(w_hat), is.numeric(S.pool), is.numeric(G.pool))

  n <- length(time)
  n0 <- sum(z == 0)
  n1 <- sum(z == 1)
  if (n0 + n1 != n) stop("z must be a (0/1) treatment group indicator")

  supported_schemes <- c("fh", "schemper", "XO", "MB", "custom_time", "fh_exp1", "fh_exp2")
  if (!(scheme %in% supported_schemes)) {
    stop("scheme must be one of: ", paste(supported_schemes, collapse = ", "))
  }
  # Only include tpoints in main list for schemes that do NOT require it in scheme_params
  if (scheme %in% c("MB", "custom_time")) {
    scheme_params$tpoints <- atpoints
    wt_args <- c(list(S = S.pool, scheme = scheme), scheme_params)
  } else {
    wt_args <- c(list(S = S.pool, scheme = scheme, tpoints = atpoints), scheme_params)
  }
  # Parameter validation for each scheme
  if (scheme == "fh" && (is.null(scheme_params$rho) || is.null(scheme_params$gamma))) {
    stop("For Fleming-Harrington weights, specify both rho and gamma in scheme_params.")
  }
  if (scheme == "schemper" && (is.null(scheme_params$Scensor) || length(scheme_params$Scensor) != length(S.pool))) {
    stop("For Schemper weights, provide Scensor (censoring KM) of same length as S in scheme_params.")
  }
  if (scheme == "XO" && (is.null(scheme_params$Ybar) || length(scheme_params$Ybar) != length(S.pool))) {
    stop("For XO weights, provide Ybar (risk set sizes) of same length as S in scheme_params.")
  }
  if (scheme == "MB" && (is.null(scheme_params$tpoints) || length(scheme_params$tpoints) != length(S.pool) || is.null(scheme_params$mb_tstar))) {
    stop("For MB weights, provide tpoints (time points) of same length as S and mb_tstar (cutoff time) in scheme_params.")
  }
  if (scheme == "custom_time" && (is.null(scheme_params$tpoints) || length(scheme_params$tpoints) != length(S.pool) || is.null(scheme_params$t.tau))) {
    stop("For custom_time weights, provide tpoints (time points) of same length as S and t.tau (cutoff time) in scheme_params.")
  }
  # Calculate weights
  wt_rg <- tryCatch(
    do.call(wt.rg.S, wt_args),
    error = function(e) {
      stop(sprintf("Weight calculation failed for scheme %s: %s", scheme, e$message))
    }
  )
  get_Cox <- tryCatch(
    uniroot(f = cox_score_rhogamma, interval = c(-15, 15), extendInt = "yes", tol = 1e-10,
            time = time, delta = delta, z = z, w_hat = w_hat, wt_rg = wt_rg),
    error = function(e) NA
  )
  if (is.na(get_Cox$root)) {
    warning("Root finding failed.")
    return(list(bhat = NA, u.beta = NA, u.zero = NA, status = "fail"))
  }
  bhat_rhogamma <- get_Cox$root
  u.zero <- cox_score_rhogamma(beta = 0, time = time, delta = delta, w_hat = w_hat, z = z, wt_rg = wt_rg)
  u.beta <- cox_score_rhogamma(beta = bhat_rhogamma, time = time, delta = delta, w_hat = w_hat, z = z, wt_rg = wt_rg)
  return(list(bhat = bhat_rhogamma, u.beta = u.beta, u.zero = u.zero, status = "ok", wt_rg = wt_rg, time = time, delta = delta, z = z, w_hat = w_hat))
}








# Required packages
#' @import survival
#' @importFrom future plan multisession
#' @importFrom future.apply future_lapply


# ---- Helper Functions ----


#' Count weighted events with delta up to time x
#' @param x Time point
#' @param error Event times
#' @param delta Event indicator
#' @param weight Weights
#' @return Weighted event count
#' @export
N_rhogamma <- function(x, error, delta, weight = 1) {
  sum(weight * delta * (error <= x))
}


#' Validate scheme parameters for weighted Cox model
validate_scheme_params <- function(scheme, scheme_params, S.pool, tpoints) {
  if (scheme == "fh" && (is.null(scheme_params$rho) || is.null(scheme_params$gamma))) {
    stop("For Fleming-Harrington weights, specify both rho and gamma in scheme_params.")
  }
  if (scheme == "schemper" && (is.null(scheme_params$Scensor) || length(scheme_params$Scensor) != length(S.pool))) {
    stop("For Schemper weights, provide Scensor (censoring KM) of same length as S in scheme_params.")
  }
  if (scheme == "XO" && (is.null(scheme_params$Ybar) || length(scheme_params$Ybar) != length(S.pool))) {
    stop("For XO weights, provide Ybar (risk set sizes) of same length as S in scheme_params.")
  }
  if (scheme == "MB" && (is.null(scheme_params$tpoints) || length(scheme_params$tpoints) != length(S.pool) || is.null(scheme_params$mb_tstar))) {
    stop("For MB weights, provide tpoints (time points) of same length as S and mb_tstar (cutoff time) in scheme_params.")
  }
  if (scheme == "custom_time" && (is.null(scheme_params$tpoints) || length(scheme_params$tpoints) != length(S.pool) || is.null(scheme_params$t.tau))) {
    stop("For custom_time weights, provide tpoints (time points) of same length as S and t.tau (cutoff time) in scheme_params.")
  }
}


# Only include tpoints in main list for schemes that do NOT require it in scheme_params
# if (scheme %in% c("MB", "custom_time")) {
#   scheme_params$tpoints <- tpoints
#   wt_args <- c(list(S = S.pool, scheme = scheme), scheme_params)
# } else {
#   wt_args <- c(list(S = S.pool, scheme = scheme, tpoints = tpoints), scheme_params)
# }



#' Calculate weights for weighted Cox model
get_weights <- function(scheme, scheme_params, S.pool, tpoints) {
  if (scheme %in% c("MB", "custom_time")) {
    scheme_params$tpoints <- tpoints
    wt_args <- c(list(S = S.pool, scheme = scheme), scheme_params)
  } else {
    wt_args <- c(list(S = S.pool, scheme = scheme, tpoints = tpoints), scheme_params)
  }
  do.call(wt.rg.S, wt_args)
}

#' Root-finding for Cox score function
find_cox_root <- function(time, delta, z, w_hat, wt_rg) {
  tryCatch(
    uniroot(f = cox_score_rhogamma, interval = c(-15, 15), extendInt = "yes", tol = 1e-10,
            time = time, delta = delta, z = z, w_hat = w_hat, wt_rg = wt_rg),
    error = function(e) NA
  )
}

ci_cox  <- function(bhat, se_bhat, alpha = 0.05, verbose = FALSE) {
  z <- qnorm(1 - alpha / 2)
  bhat_lower <- bhat - z * se_bhat
  bhat_upper <- bhat + z * se_bhat
  hr <- exp(bhat)
  lower <- exp(bhat_lower)
  upper <- exp(bhat_upper)
  if (verbose) {
    cat(sprintf("Hazard Ratio (HR): %.3f\n", hr))
    cat(sprintf("95%% CI: [%.3f, %.3f]\n", lower, upper))
  }
  result <- data.frame(
    beta = bhat,
    se_beta = se_bhat,
    hr = hr,
    lower = lower,
    upper = upper
  )
return(result)
  }


#' Weighted Cox model with (rho, gamma) weights
#' @param dfcount Data frame with columns: time, delta, z, w_hat, survP_all, survG_all
#' @param scheme Weighting scheme
#' @param scheme_params List of scheme parameters
#' @return List with estimated coefficient and diagnostics
cox_rhogamma <- function(dfcount, scheme = "fh", scheme_params = list(rho = 0, gamma = 0.5), draws = 0, parallel_resampling = FALSE, alpha = 0.05, verbose = FALSE, lr.digits = 4) {
  ans <- list()
  # Extract variables
  time   <- dfcount$time
  delta  <- dfcount$delta
  z      <- dfcount$z
  w_hat  <- dfcount$w_hat
  S.pool <- dfcount$survP_all
  #atpoints <- time

  # Input checks
  stopifnot(is.numeric(time), is.numeric(delta), is.numeric(z), is.numeric(w_hat), is.numeric(S.pool))
  n <- length(time)
  n0 <- sum(z == 0)
  n1 <- sum(z == 1)
  if (n0 + n1 != n) stop("z must be a (0/1) treatment group indicator")

  # Validate scheme and parameters
  supported_schemes <- c("fh", "schemper", "XO", "MB", "custom_time", "fh_exp1", "fh_exp2")
  if (!(scheme %in% supported_schemes)) {
    stop("scheme must be one of: ", paste(supported_schemes, collapse = ", "))
  }

  # if (scheme %in% c("MB", "custom_time")) {
  #   scheme_params$tpoints <- time
  #   wt_args <- c(list(S = S.pool, scheme = scheme), scheme_params)
  # } else {
  #   wt_args <- c(list(S = S.pool, scheme = scheme, tpoints = tpoints), scheme_params)
  # }

   if(scheme == "MB"){
   scheme_params <- list(mb_tstar = scheme_params$mb_tstar, tpoints = time)
   }

  validate_scheme_params(scheme, scheme_params, S.pool, tpoints = time)

  # Calculate weights
  wt_rg <- get_weights(scheme, scheme_params, S.pool, tpoints = time)

  #wt_rg <-  do.call(wt.rg.S, wt_args)

  # Find root of score function
  get_Cox <- find_cox_root(time, delta, z, w_hat, wt_rg)

  if (is.na(get_Cox$root)) {
    warning("Root finding failed.")
    return(list(bhat = NA, u.beta = NA, u.zero = NA, status = "fail"))
  }

  bhat_rhogamma <- get_Cox$root
  u.zero <- cox_score_rhogamma(beta = 0, time = time, delta = delta, w_hat = w_hat, z = z, wt_rg = wt_rg)
  u.beta <- cox_score_rhogamma(beta = bhat_rhogamma, time = time, delta = delta, w_hat = w_hat, z = z, wt_rg = wt_rg)

  fit_rhogamma <- list(
    bhat = bhat_rhogamma,
    u.beta = u.beta,
    u.zero = u.zero,
    status = "ok",
    wt_rg = wt_rg,
    time = time,
    delta = delta,
    z = z,
    w_hat = w_hat
  )

  ans$fit <- fit_rhogamma

  get_resamples <- cox_rhogamma_resample(fit_rhogamma = fit_rhogamma, draws = draws, parallel = parallel_resampling)

  # Score test
  z.score <- u.zero / sqrt(get_resamples$sig2U.bzero)
  ans$z.score <- z.score

  pval <- 1 - pnorm(z.score)
  ans$zlogrank_text <- paste0("logrank (1-sided) p = ", format_pval(pval, eps = 0.001, digits = lr.digits))
  if(verbose) cat("z-statistic: ", ans$zlogrank_text, "\n")
  ans$fit_resamples <- get_resamples
  se_bhat_asy <- get_resamples$se.beta.asy
  # CI based on asymptotic SE
  hr_ci_asy <- ci_cox(bhat = fit_rhogamma$bhat, se_bhat = se_bhat_asy, alpha = alpha, verbose = verbose)
  ans$hr_ci_asy <- hr_ci_asy
 if(verbose)
  if(draws > 0){
  se_bhat_star <- get_resamples$se.beta
  # De-biased hr
  bhat_debiased <- fit_rhogamma$bhat - mean(get_resamples$bhat.center.star, na.rm =TRUE)
  ans$fit$bhat_debiased <- bhat_debiased
  # CI based on asymptotic SE
  ans$hr_ci_star <- ci_cox(bhat = bhat_debiased, se_bhat = se_bhat_star, alpha = alpha, verbose = verbose)
  }
ans
}



# ---- Resampling Function ----

#' Resampling for Weighted Cox Model (rho, gamma)
#' @param bhat Estimated coefficient
#' @param time Event/censoring times
#' @param delta Event indicator
#' @param z Treatment group indicator (0/1)
#' @param w_hat Subjects' weighting (eg. propensity-scores)
#' @param G1.draws, G0.draws Optional: pre-generated random draws for groups
#' @param draws Number of resampling iterations
#' @param rho, gamma Weighting parameters
#' @param t.tau, w0.tau, w1.tau Optional: custom weights
#' @param seed.value Optional: random seed
#' @param parallel Logical: use parallelization? (default FALSE)
#' @param workers Number of parallel workers
#' @return List with resampling results
#' @export
cox_rhogamma_resample <- function(fit_rhogamma,G1.draws = NULL, G0.draws = NULL,
                                  draws = 100, seedstart=8316951,
                                  parallel = FALSE, workers = NULL
) {

  bhat <- fit_rhogamma$bhat

  time <- fit_rhogamma$time
  delta <- fit_rhogamma$delta
  z <- fit_rhogamma$z
  w_hat <- fit_rhogamma$w_hat
  wt_rg <- fit_rhogamma$wt_rg
  at_points <- time

  stopifnot(is.numeric(bhat), is.numeric(time), is.numeric(delta), is.numeric(z), is.numeric(w_hat), is.numeric(wt_rg))

  n <- length(time)
  n0 <- sum(z == 0)
  n1 <- sum(z == 1)
  if (n0 + n1 != n) stop("z must be a (0/1) treatment group indicator")

  if (is.null(G0.draws) && is.null(G1.draws) && draws > 0) {
    set.seed(seedstart)
    G0.draws <- matrix(rnorm(draws * n0), ncol = draws)
    G1.draws <- matrix(rnorm(draws * n1), ncol = draws)
  }
  idx0 <- which(z == 0)
  idx1 <- which(z == 1)
  y0 <- time[idx0]; d0 <- delta[idx0]
  #w0 <- wt_rg[idx0]
  y1 <- time[idx1]; d1 <- delta[idx1]
  #w1 <- wt_rg[idx1]
  w0_hat <- w_hat[idx0]; w1_hat <- w_hat[idx1]

  event_mat0 <- outer(y0, at_points, FUN = "<=")
  risk_mat0  <- outer(y0, at_points, FUN = ">=")
  event_mat1 <- outer(y1, at_points, FUN = "<=")
  risk_mat1  <- outer(y1, at_points, FUN = ">=")
  risk_z0 <- colSums(risk_mat0 *  w0_hat)
  risk_z1 <- colSums(risk_mat1 *  w1_hat * exp(bhat))
  counting0 <- colSums(event_mat0 * (d0 *  w0_hat))
  counting1 <- colSums(event_mat1 * (d1 *  w1_hat))
  dN_z0 <- diff(c(0, counting0))
  dN_z1 <- diff(c(0, counting1))
  num <- wt_rg * risk_z1 * risk_z0
  den <- risk_z0 + risk_z1
  K <- ifelse(den > 0, num / den, 0.0)
  drisk1 <- ifelse(risk_z1 > 0, dN_z1 / risk_z1, 0.0)
  drisk0 <- ifelse(risk_z0 > 0, dN_z0 / risk_z0, 0.0)
  score_obs <- sum(K * (drisk0 - drisk1))
  i_bhat <- sum(ifelse(den > 0, (num / (den^2)) * (dN_z0 + dN_z1), 0.0))
  DNbar <- dN_z0 + dN_z1
  h1 <- ifelse(risk_z1 > 0, (K^2 / risk_z1), 0.0)
  h2 <- ifelse(risk_z0 > 0, (K^2 / risk_z0), 0.0)
  temp <- c(den - 1)
  ybar_mod <- ifelse(temp < 1, 1, temp)
  dH1 <- ifelse(ybar_mod > 0, (DNbar-1) / ybar_mod, 0.0)
  dH2 <- ifelse(den > 0, DNbar / den, 0.0)
  sig2s <- (h1+h2)*(1-dH1)*dH2
  sig2U_bzero <- sum(sig2s)
  sig_beta_asy <- sqrt(sig2U_bzero/i_bhat^2)
  if(draws > 0){
    if(!parallel){
      counting0_star_all <- t(event_mat0 *  w0_hat) %*% (d0 * G0.draws)
      counting1_star_all <- t(event_mat1 *  w1_hat) %*% (d1 * G1.draws)
      dN_z0_star_all <- apply(counting0_star_all, 2, function(x) diff(c(0, x)))
      dN_z1_star_all <- apply(counting1_star_all, 2, function(x) diff(c(0, x)))
      drisk1_star <- sweep(dN_z1_star_all, 1, risk_z1, "/")
      drisk1_star[is.infinite(drisk1_star) | is.nan(drisk1_star)] <- 0
      drisk0_star <- sweep(dN_z0_star_all, 1, risk_z0, "/")
      drisk0_star[is.infinite(drisk0_star) | is.nan(drisk0_star)] <- 0
      score_star <- colSums(K * (drisk0_star - drisk1_star))
      bhat_center_star <- score_star / i_bhat
      Z_bstar <- bhat_center_star / sig_beta_asy
      var_bhat <- var(bhat_center_star, na.rm = TRUE)
      se_beta <- sqrt(var_bhat)
      return(list(
        se.beta.asy = sig_beta_asy,
        sig2U.bzero = sig2U_bzero,
        score.star = score_star,
        bhat.center.star = bhat_center_star,
        Z.bstar = Z_bstar,
        var.bhat = var_bhat,
        se.beta = se_beta,
        i.bhat = i_bhat,
        score.obs = score_obs
      ))
    }
    if (parallel) {
      if(draws >= 50000) stop("parallel does not currently support 50k draws")
      if (!requireNamespace("future.apply", quietly = TRUE)) stop("Please install the 'future.apply' package.")
      if (is.null(workers)) {
        workers <- max(1, parallel::detectCores() - 1)
      }
      old_plan <- future::plan()
      future::plan(future::multisession, workers = workers)
      on.exit(future::plan(old_plan), add = TRUE)
      resample_fun <- future.apply::future_lapply
    } else {
      resample_fun <- lapply
    }
    resample_results <- resample_fun(1:draws, function(dd) {
      g0 <- G0.draws[, dd]
      g1 <- G1.draws[, dd]
      counting0_star <- vapply(at_points, N_rhogamma, numeric(1), error = time[z == 0], delta = delta[z == 0] * g0 * w0_hat)
      dN_z0_star <- diff(c(0, counting0_star))
      counting1_star <- vapply(at_points, N_rhogamma, numeric(1), error = time[z == 1], delta = delta[z == 1] * g1 * w1_hat)
      dN_z1_star <- diff(c(0, counting1_star))
      drisk1_star <- ifelse(risk_z1 > 0, dN_z1_star / risk_z1, 0.0)
      drisk0_star <- ifelse(risk_z0 > 0, dN_z0_star / risk_z0, 0.0)
      score_star <- sum(K * (drisk0_star - drisk1_star))
      bhat_center_star <- score_star / i_bhat
      Z_bstar <- bhat_center_star/sig_beta_asy
      list(score.star = score_star,
           bhat.center.star = bhat_center_star,
           Z.bstar = Z_bstar)
    })
    score_star_vec <- sapply(resample_results, function(x) x$score.star)
    bhat_center_star_vec <- sapply(resample_results, function(x) x$bhat.center.star)
    Z_bstar_vec <- sapply(resample_results, function(x) x$Z.bstar)
    var_bhat <- var(bhat_center_star_vec, na.rm = TRUE)
    se_beta <- sqrt(var_bhat)
    return(list(
      se.beta.asy = sig_beta_asy,
      sig2U.bzero = sig2U_bzero,
      score.star = score_star_vec,
      bhat.center.star = bhat_center_star_vec,
      Z.bstar = Z_bstar_vec,
      var.bhat = var_bhat,
      se.beta = se_beta,
      i.bhat = i_bhat,
      score.obs = score_obs
    ))
  }
  if(draws == 0){
    return(list(
      se.beta.asy = sig_beta_asy,
      sig2U.bzero = sig2U_bzero,
      score.obs = score_obs
    ))
  }
}

