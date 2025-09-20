# This is a working version
stop("weightedCox_coxph.R functions are not ready")


# weighted Cox based on coxph (Therneau's example)
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
