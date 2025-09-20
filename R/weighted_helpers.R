#' Compute Time-Dependent Weights for Survival Analysis
#'
#' Calculates weights for use in weighted log-rank and related tests, supporting standard and custom schemes:
#' Fleming-Harrington, Schemper, Xu & O'Quigley (XO), Maggir-Burman (MB), custom time-based, and exponential variants of Fleming-Harrington weights.
#'
#' @param S Numeric vector of survival probabilities (Kaplan-Meier curve).
#' @param scheme Character string specifying weighting scheme. One of:
#'   "fh" (Fleming-Harrington), "schemper", "XO", "MB", "custom_time", "fh_exp1", or "fh_exp2".
#' @param rho Numeric weighting parameter (used for "fh").
#' @param gamma Numeric weighting parameter (used for "fh").
#' @param Scensor Optional numeric vector of censoring probabilities (for "schemper").
#' @param Ybar Optional numeric vector of risk set sizes (for "XO").
#' @param tpoints Optional numeric vector of time points (for "custom_time" and "MB").
#' @param t.tau Optional time cutoff for "custom_time" weights.
#' @param w0.tau, w1.tau Weights before/after t.tau (for "custom_time").
#' @param mb_tstar Optional time for MB weights.
#' @param details Logical; if TRUE, return detailed output.
#'
#' @return Numeric vector of weights (or list if details=TRUE).
#'
#' @details
#' - "fh": Fleming-Harrington weights, requires `rho` and `gamma`.
#' - "schemper": Schemper weights, requires `Scensor`.
#' - "XO": Xu & O'Quigley weights, requires `Ybar`.
#' - "MB": Maggir-Burman weights, requires `tpoints` and `mb_tstar`.
#' - "custom_time": Custom time-based weights, requires `tpoints` and `t.tau`.
#' - "fh_exp1": Exponential variant of FH weights.
#' - "fh_exp2": Alternative exponential variant of FH weights.
#'
#' @examples
#' # Fleming-Harrington weights
#' wt.rg.S(S = surv, scheme = "fh", rho = 1, gamma = 1)
#' # Schemper weights
#' wt.rg.S(S = surv, scheme = "schemper", Scensor = censor_surv)
#' # XO weights
#' wt.rg.S(S = surv, scheme = "XO", Ybar = risk_set)
#' # MB weights
#' wt.rg.S(S = surv, scheme = "MB", tpoints = times, mb_tstar = 12)
#' # Custom time-based weights
#' wt.rg.S(S = surv, scheme = "custom_time", tpoints = times, t.tau = 6, w0.tau = 0, w1.tau = 1)
#' # Exponential FH weights
#' wt.rg.S(S = surv, scheme = "fh_exp1")
#' wt.rg.S(S = surv, scheme = "fh_exp2")
#'
#' @export
wt.rg.S <- function(
    S,
    scheme = c("fh", "schemper", "XO", "MB", "custom_time", "fh_exp1","fh_exp2"),
    rho = NULL,
    gamma = NULL,
    Scensor = NULL,
    Ybar = NULL,
    tpoints = NULL,
    t.tau = NULL,
    w0.tau = 0,
    w1.tau = 1,
    mb_tstar = NULL,
    details = FALSE
) {
  scheme <- match.arg(scheme)
  n <- length(S)
  if (!is.numeric(S) || n < 2) stop("S must be a numeric vector of survival probabilities (length >= 2).")
  S_left <- c(1, S[-n])
  wt <- rep(1, n)

  if (scheme == "fh") {
    if (is.null(rho) || is.null(gamma)) stop("For Fleming-Harrington weights, specify both rho and gamma.")
    wt <- S_left^rho * (1 - S_left)^gamma
  } else if (scheme == "schemper") {
    if (is.null(Scensor) || length(Scensor) != n) stop("For Schemper weights, provide Scensor (censoring KM) of same length as S.")
    Scensor_left <- c(1, Scensor[-n])
    wt <- ifelse(Scensor_left > 0, S_left / Scensor_left, 0)
  } else if (scheme == "XO") {
    if (is.null(Ybar) || length(Ybar) != n) stop("For XO weights, provide Ybar (risk set sizes) of same length as S.")
    wt <- ifelse(Ybar > 0, S_left / Ybar, 0)
  } else if (scheme == "MB") {
    if (is.null(tpoints) || length(tpoints) != n) stop("For MB weights, provide tpoints (time points) of same length as S.")
    if (is.null(mb_tstar)) stop("For MB weights, provide mb_tstar (cutoff time).")
    loc_tstar <- which.max(tpoints > mb_tstar)
    Shat_tzero <- if (mb_tstar <= max(tpoints)) S_left[loc_tstar] else 0.0
    mS <- pmax(S_left, Shat_tzero)
    wt <- 1 / mS
  } else if (scheme == "custom_time") {
    if (is.null(tpoints) || length(tpoints) != n) stop("For custom_time weights, provide tpoints (time points) of same length as S.")
    if (is.null(t.tau)) stop("For custom_time weights, provide t.tau (cutoff time).")
    wt <- ifelse(tpoints <= t.tau, w0.tau, w1.tau)
  } else if (scheme == "fh_exp1") {
    wt <- exp(S_left^0.5 * (1 - S_left)^0.5)
  } else if (scheme == "fh_exp2") {
    wt05 <- exp(S_left^0.5 * (1 - S_left)^0.5)
    wmax <- max(wt05)
    wt01 <- (S_left^0) * (1 - S_left)^1
    wt <- pmin(exp(wt01), wmax)
  }
  else {
    stop("Unknown weighting scheme.")
  }
  if (details) {
    return(list(weights = wt, S = S, S_left = S_left, Scensor = Scensor, Ybar = Ybar, tpoints = tpoints, scheme = scheme))
  } else {
    return(wt)
  }
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



extract_and_calc_weights <- function(atpoints, S.pool, weights_spec_list) {
  # weights_spec_list: named list, each element is a list of scheme_params (must include 'scheme')
  # Returns: data.frame with columns: time, weight, label, scheme

  df_all <- do.call(rbind, lapply(names(weights_spec_list), function(lbl) {
    scheme_params <- weights_spec_list[[lbl]]
    if (is.null(scheme_params$scheme)) stop(paste("Missing 'scheme' in weights_spec_list for", lbl))
    scheme <- scheme_params$scheme
    scheme_params$scheme <- NULL # Remove scheme from params for get_weights
    wt <- get_weights(scheme, scheme_params, S.pool, tpoints = atpoints)
    data.frame(time = atpoints, weight = wt, label = lbl, scheme = scheme, stringsAsFactors = FALSE)
  }))
  rownames(df_all) <- NULL
  return(df_all)
}


plot_weight_schemes <- function(
    dfcount,
    tte.name = "time_months",
    event.name = "status",
    treat.name = "treat",
    arms = c("treat", "control"),
    weights_spec_list = list(
      "MB(12)"   = list(scheme = "MB", mb_tstar = 12),
      "MB(6)"    = list(scheme = "MB", mb_tstar = 6),
      "FH(0,1)"  = list(scheme = "fh", rho = 0, gamma = 1),
      "FH(0.5,0.5)"  = list(scheme = "fh", rho = 0.5, gamma = 0.5),
      "FHexp2" = list(scheme="fh_exp2")
    ),
    custom_colors = c(
      "FH(0,1)" = "grey",
      "FHexp2" = "black",
      "MB(12)" = "#1b9e77",
      "MB(6)" = "#d95f02",
      "FH(0.5,0.5)" = "#7570b3"
    ),
    custom_sizes = c(
      "FH(0,1)" = 2,
      "FHexp2" = 1,
      "MB(12)" = 1,
      "MB(6)" = 1,
      "FH(0.5,0.5)" = 1
    ),
    transform_fh = TRUE,
    rescheme_fhexp2 = FALSE,
    save_plot = FALSE,
    filename = "weights_plot.png"
) {
  # Extract atpoints, S.pool, ybar from dfcount
  atpoints <- dfcount$at_points_all
  S.pool <- dfcount$survP_all
  ybar <- dfcount$ybar_all

  # Calculate weights
  df_weights <- extract_and_calc_weights(atpoints, S.pool, weights_spec_list)


  if(!rescheme_fhexp2) df_weights$facet_group <- ifelse(df_weights$scheme == 'MB', 'MB', 'FH/FHexp2')
  if(rescheme_fhexp2) df_weights$facet_group <- ifelse(df_weights$scheme %in% c('MB','fh_exp2'), 'MB/FHexp2', 'FH')

  # Optionally transform weights for 'fh' schemes
  df_weights$weight_trans <- with(df_weights,
                                  if (transform_fh) ifelse(scheme %in% c("fh"), exp(weight), weight) else weight
  )


  # Plot
  library(ggplot2)
  g <- ggplot(df_weights, aes(x = time, y = weight_trans, linetype = label)) +
    geom_line(aes(color = label, size = label)) +
    facet_wrap(~ facet_group, scales = 'free_y') +
    scale_color_manual(values = custom_colors) +
    scale_size_manual(values = custom_sizes) +
    labs(x = 'Time (months)', y = 'Weight', linetype = 'Label',
         title = 'MB vs FH and exponential FH variant weights') +
    theme_minimal()

  if (save_plot) {
    ggsave(filename, g, width = 10, height = 6)
  }
  return(g)
}


