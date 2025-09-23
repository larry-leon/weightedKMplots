# Mean for general S (truncated at L)
# \int_{0}^{L}S(t)dt
# fatx=S, x=t
LS.int<-function(fatx,x,L=Inf){
  f <- fatx[1:(length(fatx)-1)]
  d.x <- diff(x)
  int.f <- x[1]+sum((d.x*f)[which(x[-1]<=L)])
  return(int.f)
}


mu.L<-function(L,S,t){
  mu <- LS.int(fatx=S,x=t,L=L)
  return(mu)
}

# Need to check!
LS.int.del<-function(fx,x,dx,L=Inf){
  int.f < -sum((dx*fx)[which(x<=L)])
  return(int.f)
}

# Approaches for calculating RMST (Note: implemented last version)
# Assume tpoints is sorted and dhat is the difference at each tpoint
#dt <- c(0, diff(tpoints))  # time increments
#cum_area <- cumsum(dhat * dt)  # cumulative area (trapezoidal rule can be used for more accuracy)
# dt <- diff(tpoints)
# mid_dhat <- (head(dhat, -1) + tail(dhat, -1)) / 2
# cum_area_trapz <- cumsum(mid_dhat * dt)
# cum_area_trapz <- c(0, cum_area_trapz)  # align with tpoints
# Vectorized cumulative integration for all draws
# dt <- diff(tpoints)
# mid_dhat_star <- (head(dhat_star_mat, -1) + tail(dhat_star_mat, -1)) / 2
# cum_area_stars <- apply(mid_dhat_star, 2, function(x) cumsum(x * dt))
# cum_area_stars <- rbind(0, cum_area_stars)  # align with tpoints
# cum_area_stars <- apply(dhat_star_mat, 2, function(dhat_star) {
#   mid_dhat <- (head(dhat_star, -1) + tail(dhat_star, -1)) / 2
#   c(0, cumsum(mid_dhat * dt))
# })


cumulative_rmst_bands <- function(df, fit, tte.name, event.name, treat.name, weight.name = NULL, draws_sb = 1000, xlab="months", ylim_pad = 0.5,
                                  rmst_max_legend = "left", rmst_max_cex = 0.7){
ans <- list()
at_points <- fit$at_points
dhat <- fit$dhat
dhat_star_mat <- fit$dhat_star
risk.points <- fit$risk.points

ans$at_points <- at_points

dt <- diff(tpoints)
mid_dhat <- (head(dhat, -1) + tail(dhat, -1)) / 2

rmst_time <- c(0, cumsum(mid_dhat * dt))

ans$rmst_time <- rmst_time

# For multiple draws (matrix: rows=time, cols=draws)
cum_rmst_stars <- apply(dhat_star_mat, 2, function(dhat_star) {
  mid_dhat <- (head(dhat_star, -1) + tail(dhat_star, -1)) / 2
  c(0, cumsum(mid_dhat * dt))
})
sig2_rmst_time <- apply(cum_rmst_stars, 1, var, na.rm=TRUE)

ans$sig2_rmst_time <- sig2_rmst_time

# Pointwise CIs
rmst_time_lower <- rmst_time - 1.96*sqrt(sig2_rmst_time)
rmst_time_upper <- rmst_time + 1.96*sqrt(sig2_rmst_time)

ans$rmst_time_lower <- rmst_time_lower
ans$rmst_time_upper <- rmst_time_upper

# Extract RMST at max time
at_maxtau <- length(rmst_time)
rmst_maxtau <- rmst_time[at_maxtau]
rmst_maxtau_lower <- rmst_time_lower[at_maxtau]
rmst_maxtau_upper <- rmst_time_upper[at_maxtau]
rmst_maxtau_ci <- c(rmst_maxtau, rmst_maxtau_lower, rmst_maxtau_upper)

ans$rmst_maxtau_ci <- rmst_maxtau_ci

rmst_text <- paste0("RMST(tau*) = ", round(rmst_maxtau, 1),
                   " (", round(rmst_maxtau_lower, 1), ", ", round(rmst_maxtau_upper, 1), ")")

ans$rmst_text <- rmst_text

# Centered resamples
fit_draws <- KM_diff(
  df = df, tte.name = tte.name, event.name = event.name , treat.name = treat.name, weight.name = weight.name,
  at_points = at_points, alpha = 0.05, risk.points = risk.points,
  modify_tau = FALSE,
  draws.band = draws_sb, seedstart = 99999, show_resamples = FALSE
)

dhat_star_mat2 <- fit_draws$dhat_star

cum_rmst_stars2 <- apply(dhat_star_mat2, 2, function(dhat_star) {
  mid_dhat <- (head(dhat_star, -1) + tail(dhat_star, -1)) / 2
  c(0, cumsum(mid_dhat * dt))
})

# Standardized resamples
rmst_time_star <- cum_rmst_stars2 / sqrt(sig2_rmst_time)

# simultaneous band
sups <- apply(abs(rmst_time_star), 2, max, na.rm = TRUE)
c_alpha_band <- quantile(sups,c(0.95))

ans$c_alpha_band <- c_alpha_band

# Simultaneous bands
rmst_time_sb_lower <- rmst_time - c_alpha_band * sqrt(sig2_rmst_time)
rmst_time_sb_upper <- rmst_time + c_alpha_band * sqrt(sig2_rmst_time)

ans$rmst_time_sb_lower <- rmst_time_sb_lower
ans$rmst_time_sb_upper <- rmst_time_sb_upper

x <- at_points
mean.value <- rmst_time
l0_pw <- rmst_time_lower
u0_pw <- rmst_time_upper
l0_sb <- rmst_time_sb_lower
u0_sb <- rmst_time_sb_upper

time.zero.label <- 0.0

ymin <- min(l0_sb)
ymax <- max(u0_sb) + ylim_pad

plot(
  x[order(x)], mean.value[order(x)], type = "n", xlab = xlab, lty = 1,
  ylab = "Cumulative RMST", ylim = c(ymin,ymax), cex.lab = 1
)
polygon(
  c(x[order(x)], rev(x[order(x)])),
  c(l0_sb[order(x)], rev(u0_sb[order(x)])),
  col = "lightgrey", border = FALSE
)
lines(x[order(x)], l0_pw[order(x)], lty=2, type="s")
lines(x[order(x)], u0_pw[order(x)], lty=2, type="s")

lines(x[order(x)], mean.value[order(x)], lty = 1, lwd = 1, type = "s")
abline(h = time.zero.label, lty = 1, col = "blue", lwd = 0.5)

legend(rmst_max_legend, legend = rmst_text, cex = rmst_max_cex, bty = "n")


return(invisible(ans))
}

