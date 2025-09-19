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

ci_cox  <- function(bhat, sig_bhat, alpha = 0.05, verbose = FALSE) {
  z <- qnorm(1 - alpha / 2)
  bhat_lower <- bhat - z * sig_bhat
  bhat_upper <- bhat + z * sig_bhat
  hr <- exp(bhat)
  lower <- exp(bhat_lower)
  upper <- exp(bhat_upper)
  if (verbose) {
    cat(sprintf("Hazard Ratio (HR): %.3f\n", hr))
    cat(sprintf("95%% CI: [%.3f, %.3f]\n", lower, upper))
  }
  result <- data.frame(
    beta = bhat,
    sig_bhat = sig_bhat,
    hr = hr,
    lower = lower,
    upper = upper
  )
  return(result)
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
