# Survival Analysis and KM Plotting - Refactored Version

library(survival)
# survminer is only used for 1 plot comparison
#library(survminer)

library(weightedKMplots)

# ---- Data Preparation ----

# Prepare GBSG data
df_gbsg <- gbsg
df_gbsg$time_months <- df_gbsg$rfstime / 30.4375
tte.name <- "time_months"
event.name <- "status"
treat.name <- "hormon"
arms <- c("treat", "control")

# ---- Main Analyses ----

# 1. GBSG - ITT Analysis
dfcount_gbsg <- get_dfcounting(df=df_gbsg, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms)
# check log-rank statistics
check_results(dfcount_gbsg)

# 2. Stratified by grade
# do not save, just for checking
res <- get_dfcounting(df=df_gbsg, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, strata.name="grade")
# stratified
with(res,lr_stratified^2/sig2_lr_stratified)
# non-stratified (should be same as above)
with(res,lr^2/sig2_lr)
# survdiff stratified
res$logrank_results$chisq
# non-stratified
with(res,z.score^2)
# stratified
with(res,z.score_stratified^2)

# 3. Subgroup Analyses (meno)
for (meno_val in c(0, 1)) {
  df_sg <- subset(df_gbsg, meno == meno_val)
  sg_res <- get_dfcounting(df=df_sg, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms)
  check_results(sg_res)
}

# 4. Add duplicate subject (pid==51) for testing
df_add <- subset(df_gbsg, pid == 51)
df_add$status <- 1.0
df_mod <- rbind(df_gbsg, df_add)
dfcount_mod <- get_dfcounting(df=df_mod, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms)
check_results(dfcount_mod)

# Compare with survminer
# km_fit_mod <- survfit(Surv(time_months,status) ~ hormon, data=df_mod)
# ggsurvplot(km_fit_mod,conf.int=FALSE,risk.table = TRUE, break.time.by=12, xlim=c(0,86),
#            tables.height = 0.2,tables.theme = theme_cleantable(),censor=TRUE)



# ---- Plotting ----

par(mfrow=c(1,2))
# GBSG KM plots
plot_weighted_km(dfcount=dfcount_gbsg)
# compare with survfit
plot_km(df=df_gbsg, tte.name=tte.name, event.name=event.name, treat.name=treat.name)

# For modified data
par(mfrow=c(1,2))
plot_weighted_km(dfcount_mod)
plot_km(df=df_mod, tte.name=tte.name, event.name=event.name, treat.name=treat.name)


# ---- Subgroup Band Plot ----
par(mfrow=c(1,1))
sg_cols <- c("blue", "brown", "green", "yellow")
  temp <- plotKM.band_subgroups(
    df=df_gbsg,
    sg_labels=c("er==0", "er==1", "grade==1", "meno==0"),
    sg_colors=sg_cols,
    xlabel="Months", ylabel="Survival differences",
    yseq_length=5, cex_Yaxis=0.7, risk_cex=0.7,
    tau_add=42, by.risk=6, risk_delta=0.05, risk.pad=0.03,
    tte.name=tte.name, treat.name=treat.name, event.name=event.name
  )
  legend("topleft", c("ITT", "er=0", "er=1", "grade 1", "pre-meno"),
         lwd=2, col=c("black", sg_cols), bty="n", cex=0.75)

# ---- Propensity Score Weighting (Rotterdam) ----

# Prepare Rotterdam data
gbsg_validate <- within(rotterdam, {
  rfstime <- ifelse(recur == 1, rtime, dtime)
  t_months <- rfstime / 30.4375
  time_months <- t_months
  status <- pmax(recur, death)
  ignore <- (recur == 0 & death == 1 & rtime < dtime)
  status2 <- ifelse(recur == 1 | ignore, recur, death)
  rfstime2 <- ifelse(recur == 1 | ignore, rtime, dtime)
  time_months2 <- rfstime2 / 30.4375
  grade3 <- ifelse(grade == "3", 1, 0)
  treat <- hormon
  id <- as.numeric(1:nrow(rotterdam))
  SG0 <- ifelse(er <= 0, 0, 1)
})

tte.name <- "time_months"
event.name <- "status"
treat.name <- "treat"
wt.name <- "sw.weights"

# Node positive only
df_rotterdam <- subset(gbsg_validate, nodes > 0)

# Propensity score model
fit.ps <- glm(treat ~ age + meno + size + grade3 + nodes + pgr + chemo + er, data=df_rotterdam, family="binomial")
pihat <- fit.ps$fitted
pihat.null <- glm(treat ~ 1, family="binomial", data=df_rotterdam)$fitted
wt.1 <- pihat.null / pihat
wt.0 <- (1 - pihat.null) / (1 - pihat)
df_rotterdam$sw.weights <- ifelse(df_rotterdam$treat == 1, wt.1, wt.0)

# Weighted and unweighted analyses
dfcount_rotterdam_unwtd <- get_dfcounting(df=df_rotterdam, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24)
dfcount_rotterdam_wtd <- get_dfcounting(df=df_rotterdam, tte.name=tte.name, event.name=event.name, treat.name=treat.name, arms=arms, by.risk=24, weight.name=wt.name, check.KM=FALSE)

# ---- Plotting Weighted vs Unweighted ----

par(mfrow=c(1,2))
plot_weighted_km(dfcount=dfcount_rotterdam_unwtd)
plot_weighted_km(dfcount=dfcount_rotterdam_wtd)

# Compare with GBSG trial data
par(mfrow=c(1,2))
plot_weighted_km(dfcount=dfcount_gbsg)
plot_weighted_km(dfcount=dfcount_rotterdam_wtd)

# ---- End of Script ----
