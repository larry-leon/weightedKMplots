
# score calculation
score_calculation <- function(ybar1, ybar0, dN1, dN0, wt_rg){
dN <- dN1 + dN0
num <- wt_rg * ybar1 * ybar0
den <- ybar0 + ybar1
K <- ifelse(den > 0, num / den, 0.0)

drisk1 <- ifelse(ybar1 > 0, dN1 / ybar1, 0.0)
drisk0 <- ifelse(ybar0 > 0, dN0 / ybar0, 0.0)
score <- sum(K * (drisk0 - drisk1))

h1 <- ifelse(ybar1 > 0, (K^2 / ybar1), 0.0)
h2 <- ifelse(ybar0 > 0, (K^2 / ybar0), 0.0)
temp <- c(den - 1)
ybar_mod <- ifelse(temp < 1, 1, temp)
dH1 <- ifelse(ybar_mod > 0, (dN-1) / ybar_mod, 0.0)
dH2 <- ifelse(den > 0, dN / den, 0.0)
sig2s <- (h1+h2)*(1-dH1)*dH2
sig2U <- sum(sig2s)
i_bhat <- sum(ifelse(den > 0, (num / (den^2)) * (dN0 + dN1), 0.0))
sig2_beta_asy <- sig2U/i_bhat^2
return(list(score = score, sig2U = sig2U, sig2_beta_asy = sig2_beta_asy, i_bhat = i_bhat, K_wt_rg = K))
}


cox_score_rhogamma <- function(beta, time, delta, z, w_hat = rep(1,length(time)), wt_rg = rep(1,length(time)), score_only = TRUE) {
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

  score_stats <- score_calculation(ybar1 = risk_z1, ybar0 = risk_z0, dN0 = dN_z0, dN1 = dN_z1, wt_rg = wt_rg)

  score <- score_stats$score

  if(score_only){
    return(score)
  }
  else{
  sig2U <- score_stats$sig2U
  sig2_beta_asy <- score_stats$sig2_beta_asy
  K_wt_rg <- score_stats$K_wt_rg
  i_bhat <- score_stats$i_bhat
  return(list(score = score, sig2_score = sig2U, sig2_beta_asy = sig2_beta_asy, K_wt_rg = K_wt_rg, i_bhat = i_bhat))
  }
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

   if(scheme == "MB"){
  if(is.null(scheme_params$mb_tstar)){
   cat("Missing mb_tstar argument in scheme_params you have:", paste(names(scheme_params), collapse = ", "), "\n")
     }
   scheme_params <- list(mb_tstar = scheme_params$mb_tstar, tpoints = time)
   }

  validate_scheme_params(scheme, scheme_params, S.pool, tpoints = time)

  # Calculate weights
  wt_rg <- get_weights(scheme, scheme_params, S.pool, tpoints = time)

  # Find root of score function
  get_Cox <- find_cox_root(time, delta, z, w_hat, wt_rg)

  if (is.na(get_Cox$root)) {
    warning("Root finding failed.")
    return(list(bhat = NA, u.beta = NA, u.zero = NA, status = "fail"))
  }

  bhat_rhogamma <- get_Cox$root

  # Score test: U(beta=0)
  temp <- cox_score_rhogamma(beta = 0, time = time, delta = delta, w_hat = w_hat, z = z, wt_rg = wt_rg, score_only = FALSE)
  u.zero <- temp$score
  z.score <- u.zero / sqrt(temp$sig2_score)
  i_zero <- temp$i_bhat
  K_zero <- temp$K_wt_rg
  sig2_score <- temp$sig2_score
  rm("temp")

  ans$z.score <- z.score

  temp <- cox_score_rhogamma(beta = bhat_rhogamma, time = time, delta = delta, w_hat = w_hat, z = z, wt_rg = wt_rg, score_only = FALSE)
  u.beta <- temp$score
  sig_bhat_asy <- sqrt(temp$sig2_beta_asy)
  i_bhat <- temp$i_bhat
  K_wt_rg <- temp$K_wt_rg
  rm("temp")

  pval <- 1 - pnorm(z.score)
  ans$zlogrank_text <- paste0("logrank (1-sided) p = ", format_pval(pval, eps = 0.001, digits = lr.digits))
  if(verbose) cat("z-statistic: ", ans$zlogrank_text, "\n")

  # CI based on asymptotic SE
  hr_ci_asy <- ci_cox(bhat = bhat_rhogamma, sig_bhat = sig_bhat_asy, alpha = alpha, verbose = verbose)
  ans$hr_ci_asy <- hr_ci_asy

  fit_rhogamma <- list(
    bhat = bhat_rhogamma,
    sig_bhat_asy = sig_bhat_asy,
    u.beta = u.beta,
    u.zero = u.zero,
    z.score = z.score,
    sig2_score = sig2_score,
    status = "ok",
    wt_rg = wt_rg,
    time = time,
    delta = delta,
    z = z,
    w_hat = w_hat
  )

  ans$fit <- fit_rhogamma

  if(draws > 0){

  get_resamples <- cox_rhogamma_resample(fit_rhogamma = fit_rhogamma, i_bhat = i_bhat, K_wt_rg = K_wt_rg,
                                         i_zero = i_zero, K_zero = K_zero, draws = draws, parallel = parallel_resampling
                                         )
  ans$fit_resamples <- get_resamples
  sig_bhat_star <- get_resamples$sig_bhat_star
  ans$fit$sig_bhat_star <- sig_bhat_star
  # De-biased hr
  bhat_debiased <- fit_rhogamma$bhat - mean(get_resamples$bhat_center_star, na.rm = TRUE)
  ans$fit$bhat_debiased <- bhat_debiased
  ans$fit$wald_debiased1 <- bhat_debiased / sig_bhat_asy
  ans$fit$wald_debiased2 <- bhat_debiased / sig_bhat_star
  ans$fit$wald <- fit_rhogamma$bhat / sig_bhat_asy
  # CI based on de-biased and resampled SEs
  ans$hr_ci_star <- ci_cox(bhat = bhat_debiased, sig_bhat = sig_bhat_star, alpha = alpha, verbose = verbose)
  # De-biased score
  sig2_score_star <- var(get_resamples$score_star, na.rm = TRUE)
  ans$sig2_score_star <- sig2_score_star
  u.zero_debiased <- u.zero - mean(get_resamples$score_star_null, na.rm = TRUE)
  ans$z.score_debiased <- u.zero_debiased / sqrt(sig2_score)
  #ans$z.score_debiased <- u.zero / sqrt(sig2_score_star)

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
cox_rhogamma_resample <- function(fit_rhogamma, i_bhat, K_wt_rg, i_zero, K_zero, G1.draws = NULL, G0.draws = NULL,
                                  draws = 100, seedstart=8316951,
                                  parallel = FALSE, workers = NULL
) {

  bhat <- fit_rhogamma$bhat
  sig_bhat_asy <- fit_rhogamma$sig_bhat_asy
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
  y1 <- time[idx1]; d1 <- delta[idx1]
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

      score_star <- colSums(K_wt_rg * (drisk0_star - drisk1_star))

      bhat_center_star <- score_star / i_bhat
      Z_bstar <- bhat_center_star / sig_bhat_asy
      var_bhat_star <- var(bhat_center_star, na.rm = TRUE)
      sig_bhat_star <- sqrt(var_bhat_star)

       # Logrank stat (replace risk_z1(bhat) with risk_z1_null = risk_z1(bhat=0))
       risk_z1_null <- colSums(risk_mat1 *  w1_hat)
       drisk1_star_null <- sweep(dN_z1_star_all, 1, risk_z1_null, "/")
       drisk1_star_null[is.infinite(drisk1_star_null) | is.nan(drisk1_star_null)] <- 0
       score_star_null <- colSums(K_zero * (drisk0_star - drisk1_star_null))

      return(list(
        score_star = score_star,
        score_star_null = score_star_null,
        bhat_center_star = bhat_center_star,
        Z_bstar = Z_bstar,
        sig_bhat_star = sig_bhat_star
      ))
    }
    if (parallel) {
      stop("This needs to be re-written")
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
      sig_bhat_asy = sig_bhat_asy,
      score_obs = score_obs
    ))
  }
}

