rm(list=ls())
library(survival)
library(weightedKMplots)

# load these prior
source("R/kmplotting_helpers.R")
source("R/df_counting_improved.R")
source("R/km_calculations_helpers.R")
source("R/weightedCox.R")
source("R/weightedCox_simulations_helpers.R")


library(gsDesign2)
library(dplyr)
library(tibble)
library(gt)
library(simtrial)
library(tidyr)
library(doFuture)
library(foreach)
library(tictoc)
library(data.table)

library(future.callr)


sims_setup <- function(surv_24 = 0.35, hr_true = log(0.35) / log(0.25), control_median = 12){
survival_at_24_months <- surv_24
hr <- hr_true
control_rate <- c(log(2) / control_median, (log(.25) - log(.2)) / 12)

ans <- list()

scenarios <- tribble(
  ~Scenario, ~Name,           ~Period, ~duration, ~Survival,
  0,         "Control",       0,       0,         1,
  0,         "Control",       1,       24,        .25,
  0,         "Control",       2,       12,        .2,
  1,         "PH",            0,       0,         1,
  1,         "PH",            1,       24,        .35,
  1,         "PH",            2,       12,        .2^hr,
  2,         "3-month delay", 0,       0,         1,
  2,         "3-month delay", 1,       3,         exp(-3 * control_rate[1]),
  2,         "3-month delay", 2,       21,        .35,
  2,         "3-month delay", 3,       12,        .2^hr,
  3,         "6-month delay", 0,       0,         1,
  3,         "6-month delay", 1,       6,         exp(-6 * control_rate[1]),
  3,         "6-month delay", 2,       18,        .35,
  3,         "6-month delay", 3,       12,        .2^hr,
  4,         "Crossing",      0,       0,         1,
  4,         "Crossing",      1,       3,         exp(-3 * control_rate[1] * 1.3),
  4,         "Crossing",      2,       21,        .35,
  4,         "Crossing",      3,       12,        .2^hr,
  5,         "Weak null",     0,       0,         1,
  5,         "Weak null",     1,       24,        .25,
  5,         "Weak null",     2,       12,        .2,
  6,         "Strong null",   0,       0,         1,
  6,         "Strong null",   1,       3,         exp(-3 * control_rate[1] * 1.5),
  6,         "Strong null",   2,       3,         exp(-6 * control_rate[1]),
  6,         "Strong null",   3,       18,        .25,
  6,         "Strong null",   4,       12,        .2,
)
# scenarios |> gt()

ans$scenarios <- scenarios

fr <-
  scenarios |>
  group_by(Scenario) |>
  #  filter(Scenario == 2) |>
  mutate(Month = cumsum(duration),
         x_rate = -(log(Survival) - log(lag(Survival, default = 1))) /
           duration,
         rate = ifelse(Month > 24, control_rate[2], control_rate[1]),
         hr = x_rate / rate) |>
  select(-x_rate) |>
  filter(Period > 0, Scenario > 0) |> ungroup()

fr <- fr |> mutate(fail_rate = rate, dropout_rate =0.001, stratum = "All")

ans$fr <- fr

# MWLR
mwlr <- fixed_design_mb(
  tau = 12,
  enroll_rate = define_enroll_rate(duration = 12, rate = 1),
  fail_rate = fr |> filter(Scenario == 2),
  alpha = 0.025, power = .85, ratio = 1,
  study_duration = 36
) |> to_integer()


ans$er <- mwlr$enroll_rate


# Constant dropout rate for both treatment arms and all scenarios
dropout_rate <- data.frame(stratum = rep("All", 2),
                           period = rep(1, 2),
                           treatment = c("control", "experimental"),
                           duration = rep(100, 2),
                           rate = rep(.001, 2)
)

ans$dropout_rate <- dropout_rate
return(ans)
}

# Last n_sim was 10000:  For seedstart to append to last n_sims
n_last <- 10*1000
seedstart <- 8316951 + 1000*n_last
these_sims <- get_sims(n_sim = 20*1000, seedstart = seedstart, num_workers = 100, mart_draws = 250, save_results = TRUE, file_togo = "vignettes/results/wCoxsims1_Linux_20k.RData")

