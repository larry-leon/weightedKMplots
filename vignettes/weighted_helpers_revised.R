# Compute Time-Dependent Weights for Survival Analysis
wt.rg.S <- function(
    S,
    scheme = c('fh', 'schemper', 'XO', 'MB', 'custom_time', 'fh_exp1','fh_exp2'),
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
  if (!is.numeric(S) || n < 2) stop('S must be a numeric vector of survival probabilities (length >= 2).')
  S_left <- c(1, S[-n])
  wt <- rep(1, n)
  if (scheme == 'fh') {
    if (is.null(rho) || is.null(gamma)) stop('For Fleming-Harrington weights, specify both rho and gamma.')
    wt <- S_left^rho * (1 - S_left)^gamma
  } else if (scheme == 'schemper') {
    if (is.null(Scensor) || length(Scensor) != n) stop('For Schemper weights, provide Scensor (censoring KM) of same length as S.')
    Scensor_left <- c(1, Scensor[-n])
    wt <- ifelse(Scensor_left > 0, S_left / Scensor_left, 0)
  } else if (scheme == 'XO') {
    if (is.null(Ybar) || length(Ybar) != n) stop('For XO weights, provide Ybar (risk set sizes) of same length as S.')
    wt <- ifelse(Ybar > 0, S_left / Ybar, 0)
  } else if (scheme == 'MB') {
    if (is.null(tpoints) || length(tpoints) != n) stop('For MB weights, provide tpoints (time points) of same length as S.')
    if (is.null(mb_tstar)) stop('For MB weights, provide mb_tstar (cutoff time).')
    loc_tstar <- which.max(tpoints > mb_tstar)
    Shat_tzero <- if (mb_tstar <= max(tpoints)) S_left[loc_tstar] else 0.0
    mS <- pmax(S_left, Shat_tzero)
    wt <- 1 / mS
  } else if (scheme == 'custom_time') {
    if (is.null(tpoints) || length(tpoints) != n) stop('For custom_time weights, provide tpoints (time points) of same length as S.')
    if (is.null(t.tau)) stop('For custom_time weights, provide t.tau (cutoff time).')
    if (is.null(w0.tau)) stop('For custom_time weights, provide w0.tau (weight before t.tau).')
    if (is.null(w1.tau)) stop('For custom_time weights, provide w1.tau (weight after t.tau).')
    wt <- ifelse(tpoints <= t.tau, w0.tau, w1.tau)
  } else if (scheme == 'fh_exp1') {
    wt <- exp(S_left^0.5 * (1 - S_left)^0.5)
  } else if (scheme == 'fh_exp2') {
    wt05 <- exp(S_left^0.5 * (1 - S_left)^0.5)
    wmax <- max(wt05)
    wt01 <- (S_left^0) * (1 - S_left)^1
    wt <- pmin(exp(wt01), wmax)
  }
  else {
    stop('Unknown weighting scheme.')
  }
  if (details) {
    return(list(weights = wt, S = S, S_left = S_left, Scensor = Scensor, Ybar = Ybar, tpoints = tpoints, scheme = scheme))
  } else {
    return(wt)
  }
}

validate_scheme_params <- function(scheme, scheme_params, S.pool, tpoints) {
  if (scheme == 'fh' && (is.null(scheme_params$rho) || is.null(scheme_params$gamma))) {
    stop('For Fleming-Harrington weights, specify both rho and gamma in scheme_params.')
  }
  if (scheme == 'schemper' && (is.null(scheme_params$Scensor) || length(scheme_params$Scensor) != length(S.pool))) {
    stop('For Schemper weights, provide Scensor (censoring KM) of same length as S in scheme_params.')
  }
  if (scheme == 'XO' && (is.null(scheme_params$Ybar) || length(scheme_params$Ybar) != length(S.pool))) {
    stop('For XO weights, provide Ybar (risk set sizes) of same length as S in scheme_params.')
  }
  if (scheme == 'MB' && (is.null(scheme_params$tpoints) || length(scheme_params$tpoints) != length(S.pool) || is.null(scheme_params$mb_tstar))) {
    stop('For MB weights, provide tpoints (time points) of same length as S and mb_tstar (cutoff time) in scheme_params.')
  }
  if (scheme == 'custom_time' && (is.null(scheme_params$t.tau) || is.null(scheme_params$w0.tau) || is.null(scheme_params$w1.tau))) {
    stop('For custom_time weights, provide tpoints (time points), t.tau (cutoff time), w0.tau, and w1.tau in scheme_params.')
  }
}

get_weights <- function(scheme, scheme_params, S.pool, tpoints) {
  if (scheme %in% c('MB', 'custom_time')) {
    scheme_params$tpoints <- tpoints
    wt_args <- c(list(S = S.pool, scheme = scheme), scheme_params)
  } else {
    wt_args <- c(list(S = S.pool, scheme = scheme, tpoints = tpoints), scheme_params)
  }
  do.call(wt.rg.S, wt_args)
}

extract_and_calc_weights <- function(atpoints, S.pool, weights_spec_list) {
  df_all <- do.call(rbind, lapply(names(weights_spec_list), function(lbl) {
    scheme_params <- weights_spec_list[[lbl]]
    if (is.null(scheme_params$scheme)) stop(paste('Missing \'scheme\' in weights_spec_list for', lbl))
    scheme <- scheme_params$scheme
    scheme_params$scheme <- NULL
    wt <- get_weights(scheme, scheme_params, S.pool, tpoints = atpoints)
    data.frame(time = atpoints, weight = wt, label = lbl, scheme = scheme, stringsAsFactors = FALSE)
  }))
  rownames(df_all) <- NULL
  return(df_all)
}

plot_weight_schemes <- function(
    dfcount,
    tte.name = 'time_months',
    event.name = 'status',
    treat.name = 'treat',
    arms = c('treat', 'control'),
    weights_spec_list = list(
      'MB(12)'   = list(scheme = 'MB', mb_tstar = 12),
      'MB(6)'    = list(scheme = 'MB', mb_tstar = 6),
      'FH(0,1)'  = list(scheme = 'fh', rho = 0, gamma = 1),
      'FH(0.5,0.5)'  = list(scheme = 'fh', rho = 0.5, gamma = 0.5),
      'custom_time' = list(scheme = 'custom_time', t.tau = 25, w0.tau = 1.0, w1.tau = 1.5),
      'FHexp2' = list(scheme='fh_exp2')
    ),
    custom_colors = c(
      'FH(0,1)' = 'grey',
      'FHexp2' = 'black',
      'MB(12)' = '#1b9e77',
      'MB(6)' = '#d95f02',
      'FH(0.5,0.5)' = '#7570b3',
      'custom_time' = "green"
    ),
    custom_sizes = c(
      'FH(0,1)' = 2,
      'FHexp2' = 1,
      'MB(12)' = 1,
      'MB(6)' = 1,
      'FH(0.5,0.5)' = 1,
      'custom_time' = 1
    ),
    transform_fh = TRUE,
    rescheme_fhexp2 = FALSE,
    save_plot = FALSE,
    filename = 'weights_plot.png'
) {
  atpoints <- dfcount$at_points_all
  S.pool <- dfcount$survP_all
  ybar <- dfcount$ybar_all
  df_weights <- extract_and_calc_weights(atpoints, S.pool, weights_spec_list)
  if(!rescheme_fhexp2) df_weights$facet_group <- ifelse(df_weights$scheme == 'MB', 'MB', 'FH/FHexp2')
  if(rescheme_fhexp2) df_weights$facet_group <- ifelse(df_weights$scheme %in% c('MB','fh_exp2','custom_time'), 'MB/FHexp2/custom', 'FH')
  df_weights$weight_trans <- with(df_weights,
                                  if (transform_fh) ifelse(scheme %in% c('fh'), exp(weight), weight) else weight
  )
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