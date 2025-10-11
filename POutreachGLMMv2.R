
## Logistic regression with random intercept (OR_site) and hi*anc4 interaction
## - Recreates a clean analysis frame (dat2)
## - Bins age with oldest group = 40+
## - Sets reference categories as requested:
##     Age_bin reference = "25–29"
##     hi reference      = "3" (None)
##     anc4 reference    = "2" (<4 ANC visits)
## - Standardizes continuous covariates to aid convergence

rm(list = ls())

## Packages
library(lme4)

## 1) Load data ---------------------------------------------------------------
dat1 <- read.csv("C:/Users/salem/Dropbox/Rose_Outreach/PregOutreach.csv",
                 header = TRUE, sep = ",")

pre <- read.csv("C:/Users/salem/Dropbox/Rose_Outreach/egadob_HB_DP.csv",
                 header = TRUE, sep = ",")



## Sanity check: 1-1 merge expectation
### will throw a warning message in the merge step if something is off
d1_dup <- dat1$ANC_code[duplicated(dat1$ANC_code)]
pr_dup <- pre$ANC_code[duplicated(pre$ANC_code)]
if (length(d1_dup)) warning("dat1 has duplicated ANC_code (examples): ", paste(head(unique(d1_dup)), collapse=", "))
if (length(pr_dup)) warning("pre has duplicated ANC_code (examples): ", paste(head(unique(pr_dup)), collapse=", "))

## Merge datasets by ANC_code
POut <- merge(dat1, pre, by = "ANC_code")

## 2) Dates & derived time ----------------------------------------------------
POut$EnterDate <- as.Date(POut$EnterDate, format = "%m/%d/%Y")
POut$EnterDays <- as.numeric(POut$EnterDate - min(POut$EnterDate, na.rm = TRUE))

## 3) Analysis frame (mirror prior work) --------------------------------------
vars_keep <- c(
  "home_birthYN",
  "yrproject","yrdfc","hx_home_sum","country_plusOR","Age", 
  "MEthnicCAT","HowLongYear","Gravida","CanRead","HouseholdHunger",
  "Caesarean","OR_site","EnterDate","EnterDays","hi","anc4","NOC_MX_max",
			   "enrolEGA", "EGA_birth", "prematureY_N")


dat2 <- POut[, vars_keep]

## 4) Factors & coding --------------------------------------------------------
## Outcome as factor (No/Yes) as in prior scripts; lme4 handles 2-level factors
dat2$home_birthYN <- factor(dat2$home_birthYN, levels = c(0,1), labels = c("No","Yes"))

dat2$MEthnicCAT   <- factor(dat2$MEthnicCAT)
dat2$OR_site      <- factor(dat2$OR_site)
dat2$Caesarean    <- factor(dat2$Caesarean)
dat2$yrproject    <- factor(dat2$yrproject)
dat2$hi           <- factor(dat2$hi)    # expected levels: "3","1","2"
dat2$anc4         <- factor(dat2$anc4)  # expected levels: "1" (>=4), "2" (<4)
dat2$CanRead <- factor(dat2$CanRead, levels = c(0,1), labels = c("No","Yes"))

dat2$prematureY_N         <- factor(dat2$prematureY_N)  # expected levels: "Y" (yes, preterm_, "N" not preterm

## Collapse hi into two categories: insur vs no_insur
## - "3" = None → no_insur
## - "1" (MFUND) + "2" (Thai HI) → insur
dat2$hi2 <- ifelse(dat2$hi %in% c("1","2"), "insur", "no_insur")
dat2$hi2 <- factor(dat2$hi2, levels = c("no_insur","insur"))  # set reference



## Age bins (5-year bins with oldest collapsed to 40+)
dat2$Age_bin <- cut(dat2$Age,
                    breaks = c(15, 20, 25, 30, 35, 40, Inf),
                    right  = FALSE,
                    labels = c("15–19","20–24","25–29","30–34","35–39","40+"))

## Set reference/comparator levels
dat2$Age_bin <- relevel(dat2$Age_bin, ref = "25–29")

dat2$hi   <- factor(dat2$hi,   levels = c("3","1","2"))  # 3=None (reference), 1=MFUND, 2=Thai HI
dat2$anc4 <- factor(dat2$anc4, levels = c("2","1"))      # 2=<4 (reference), 1=≥4

## 5) Standardize continuous covariates (mean 0, SD 1) -----------------------
## create a function for centering a variable on its mean and standardizing by its sd
safe_scale <- function(x) {
  if (all(is.na(x))) return(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(x*0)          # all zeros if no spread
  as.numeric((x - m)/s)
}

## call the variables we want to standardize
cont_vars <- c("Gravida","hx_home_sum","HowLongYear","EnterDays","EGA_birth","NOC_MX_max")
## call the function to standardize those variables
dat2[cont_vars] <- lapply(dat2[cont_vars], safe_scale)


## 5.1) Preserve raw counts BEFORE any gating --------------------------------
### in this case, 'gating' is referring to the process through which some variables are applied only if
###  certain conditions are met
###  For example: one can only have had a previous C-section, if one has a gravidity > 0
dat2$Gravida_raw     <- POut$Gravida
dat2$hx_home_sum_raw <- POut$hx_home_sum



## 6A) GATED (single-sample) model -------------------------------------------

## Decompose gravidity and gate history variables to multigravids

## Indicators / decompositions (NA-safe)
##indicator for primigravida
dat2$G1 <- ifelse(is.na(dat2$Gravida_raw), NA_integer_, as.integer(dat2$Gravida_raw == 1L))
## indicator for multigravida
dat2$M  <- ifelse(is.na(dat2$Gravida_raw), NA_integer_, as.integer(dat2$Gravida_raw >= 2L))
## count beyond first pregnancy
dat2$Gplus <- ifelse(is.na(dat2$Gravida_raw), NA_real_, pmax(dat2$Gravida_raw - 1, 0))

# Robust binary for prior C-section
## looking for variations on the code (misspellings, etc.)
dat2$Caesarean01 <- as.integer(as.character(dat2$Caesarean) %in% c("Yes","Ye","1","Y","TRUE","True"))

## Gate variables: 0 for primigravids; keep NA among multigravids if missing (later will multiply by this, so if 1 = same, if 0 = no outcome - gated)
## the logic here is: if a woman is multigravid (M==1) then retain C-section or home birth history
###  if primigravid (M==0) then force this to 0 
dat2$Csec_M   <- ifelse(dat2$M == 1L, dat2$Caesarean01, 0L)
dat2$hxhome_M <- ifelse(dat2$M == 1L, dat2$hx_home_sum_raw, 0)

## Multigravid-only standardization (zeros elsewhere)
dat2$Gplus_z <- 0
dat2$hxhome_M_z <- 0

#centering on the mean and standardizing by sdcor
## but only if there are >2 valid observations
## computes mean and sd only among multigravids
## replaces standardized values back into these rows
## keeps 0s for everyone else
####   if we don't do this, then the folks who aren't multigravids are contributing to the standardized values
idx_M <- which(dat2$M == 1L & !is.na(dat2$Gplus))
if (length(idx_M) > 2) {
  mu_gp <- mean(dat2$Gplus[idx_M], na.rm = TRUE)
  sd_gp <-  sd(dat2$Gplus[idx_M],  na.rm = TRUE)
  if (is.finite(sd_gp) && sd_gp > 0) {
    dat2$Gplus_z[idx_M] <- (dat2$Gplus[idx_M] - mu_gp) / sd_gp
  }
}

##same as idx_M above, but with history of home birth here instead
idx_H <- which(dat2$M == 1L & !is.na(dat2$hxhome_M))
if (length(idx_H) > 2) {
  mu_hx <- mean(dat2$hxhome_M[idx_H], na.rm = TRUE)
  sd_hx <-  sd(dat2$hxhome_M[idx_H],  na.rm = TRUE)
  if (is.finite(sd_hx) && sd_hx > 0) {
    dat2$hxhome_M_z[idx_H] <- (dat2$hxhome_M[idx_H] - mu_hx) / sd_hx
  }
}


##### this is the original model with anc4 (anc visits categorized rather than continuous) 
#logit_mod_gated <- glmer(
#  home_birthYN ~ hi2+anc4 + Age_bin + yrproject + MEthnicCAT + CanRead + #HouseholdHunger +
#    G1 + Gplus_z + Csec_M + hxhome_M_z + HowLongYear + EnterDays + EGA_birth +
#    (1 | OR_site),
#  data   = dat2,
#  family = binomial(link = "logit"),
#  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
#)

## Fit the gated model (replace raw Gravida/hx_home_sum with gated pieces)
logit_mod_gated <- glmer(
  home_birthYN ~ hi2+ Age_bin + yrproject + MEthnicCAT + CanRead + #HouseholdHunger +
    G1 + Gplus_z + Csec_M + hxhome_M_z + HowLongYear + EnterDays + EGA_birth + NOC_MX_max +
    (1 | OR_site),
  data   = dat2,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)



#########################################################################################################
### analysis of OR sites/clusters
### calculating both the MOR (median odds ratio) and ICC, as well as CIs for them
#########################################################################################################
## -- parametric bootstrap (serial; robust on all OS) --------------------------
# Simulate new outcomes from the fitted GLMM, refit each replicate, and compute VAR/MOR/ICC.
# Use a large B (e.g., 1000–2000) for final results.
set.seed(20250908)

#beware that this will make the process run A VERY LONG TIME
B <- 2000

stat_fun <- function(fit) {
  s2 <- get_sigma2(fit)
  c(VAR = s2, MOR = mor_from_var(s2), ICC = icc_from_var(s2))
}

library(tictoc)
tic()
bs <- bootMer(
  x        = logit_mod_gated,
  FUN      = stat_fun,
  nsim     = B,
  use.u    = FALSE,       # parametric: simulate new random effects each replicate
  type     = "parametric",
  parallel = "no",        # avoid SNOW issues; change to "multicore"/"snow" only if confirmed working
  re.form  = NULL,
  verbose  = 0
)
toc()

## -- summarize ----------------------------------------------------------------
SEs    <- apply(bs$t, 2, sd, na.rm = TRUE)
VAR_CI <- pct(bs$t[, "VAR"])
MOR_CI <- pct(bs$t[, "MOR"])
ICC_CI <- pct(bs$t[, "ICC"])

mor_icc_out <- list(
  # point estimates from original fit
  sigma2 = unname(sigma2_hat),
  MOR    = unname(MOR_hat),
  ICC    = unname(ICC_hat),

  # bootstrap SEs
  SE     = setNames(SEs, names(SEs)),

  # 95% percentile CIs from bootstrap
  VAR_CI = setNames(VAR_CI, c("lower","upper")),
  MOR_CI = setNames(MOR_CI, c("lower","upper")),
  ICC_CI = setNames(ICC_CI, c("lower","upper")),

  # diagnostics
  boot_reps      = B,
  boot_converged = sum(complete.cases(bs$t[, "VAR"]))
)

print(mor_icc_out)








## 6B) STRATIFIED models ------------------------------------------------------

#there are not enough older age groups to run the exact same age bins in these models
##    Create 25+ age groups in each subset -----------------------------------
mk_age_25plus <- function(age_vec) {
  cut(age_vec,
      breaks = c(15, 20, 25, Inf),
      right  = FALSE,
      labels = c("15–19","20–24","25+")) |>
    factor(levels = c("25+","15–19","20–24"))  # set reference first
}



## Primigravid only (G == 1): no prior C-sec/home-birth terms (not defined)
dat_prim  <- droplevels(subset(dat2, Gravida_raw == 1))

dat_prim$Age_bin25p  <- mk_age_25plus(dat_prim$Age)

logit_mod_prim <- glmer(
  home_birthYN ~ hi2+ Age_bin25p + yrproject + MEthnicCAT + CanRead +
    HowLongYear + EnterDays + EGA_birth + NOC_MX_max +
    (1 | OR_site),
  data   = dat_prim,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)


## Multigravid only (G >= 2): include prior C-sec, prior home-births, and dose
dat_multi <- droplevels(subset(dat2, Gravida_raw >= 2))

dat_multi$Age_bin25p <- mk_age_25plus(dat_multi$Age)

## You can use Gplus_z (dose beyond the first). If you prefer the original standardized "Gravida",
## replace Gplus_z with Gravida in the formula (interpretation changes to per 1 SD in total gravidity).
logit_mod_multi <- glmer(
  home_birthYN ~ hi2 + Age_bin25p + yrproject + MEthnicCAT + CanRead +
    Gplus_z + Caesarean01 + hxhome_M_z + HowLongYear + EnterDays + EGA_birth + NOC_MX_max +
    (1 | OR_site),
  data   = dat_multi,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)


## (Optional) Quick sanity checks on domains
stopifnot(all(dat_prim$Caesarean01 == 0, na.rm=TRUE))
stopifnot(all(dat_prim$hx_home_sum_raw %in% c(0, NA)))











## =============================================================================
## Unified RESULTS TABLE
## -----------------------------------------------------------------------------
## What this script produces (one row per category/variable):
##   - Descriptive counts: Yes, Total, % Homebirth (for categorical predictors)
##   - CRUDE odds ratios (UORs) from the data only
##       * Categorical: 2×2 table OR vs reference level (no modeling)
##       * Continuous: per 1 SD from univariate glm() with ONLY that predictor
##   - MODEL-ADJUSTED odds ratios (AORs) from your GLMM (with (1|OR_site) + covariates)
##   - Per-unit back-translations for continuous predictors (both UOR & AOR)
##
## Assumptions:
##   - You already have in memory:
##       dat2              : analysis frame used for the overall model
##       POut              : raw data with unstandardized continuous variables
##       logit_mod_gated   : fitted GLMM with (1 | OR_site) and your covariates
##   - Outcome is a factor with levels c("No","Yes") named "home_birthYN"
##
## Edit the 'cats' and 'conts' lists to control which variables appear and units.
## =============================================================================

suppressPackageStartupMessages({
  library(lme4)
})

## ---------- Utilities ----------
esc_re  <- function(x) gsub("([][{}()+*^$.|?\\\\])", "\\\\\\1", x)
coef_df <- function(fit) as.data.frame(summary(fit)$coefficients, stringsAsFactors = FALSE)

fmt_or <- function(beta, se) {
  OR <- exp(beta); L <- exp(beta - 1.96*se); U <- exp(beta + 1.96*se)
  sprintf("%.2f (%.2f–%.2f)", OR, L, U)
}

## Counts/percents for categorical variables
level_counts <- function(data, var, outcome_var = "home_birthYN", yes_label = "Yes") {
  x <- droplevels(data[[var]])
  if (!is.factor(x)) x <- factor(x)

  tab_all <- as.data.frame(table(x), stringsAsFactors = FALSE)
  tab_yes <- as.data.frame(table(x, data[[outcome_var]]), stringsAsFactors = FALSE)
  names(tab_yes) <- c("x", "y", "Freq")
  yes <- tab_yes[tab_yes$y == yes_label, c("x", "Freq")]
  names(yes) <- c("Level", "Yes")

  out <- merge(setNames(tab_all, c("Level", "Total")), yes, by = "Level", all.x = TRUE)
  out$Yes[is.na(out$Yes)] <- 0
  out$Percent <- ifelse(out$Total > 0, 100 * out$Yes / out$Total, NA_real_)
  out
}

## Crude 2×2 OR (categorical level vs reference) with Haldane–Anscombe correction if needed
## Inputs: counts list(list(Yes=..., No=...), list(Yes=..., No=...)) for (L, Ref)
or_2x2 <- function(yL, nL, yR, nR) {
  a <- yL; b <- nL - yL; c <- yR; d <- nR - yR
  ## Haldane–Anscombe if any cell is 0
  if (any(c(a,b,c,d) == 0)) { a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5 }
  OR  <- (a/b)/(c/d)
  seL <- sqrt(1/a + 1/b + 1/c + 1/d)
  L   <- exp(log(OR) - 1.96*seL)
  U   <- exp(log(OR) + 1.96*seL)
  sprintf("%.2f (%.2f–%.2f)", OR, L, U)
}

## Univariate GLM for a standardized continuous raw variable (crude per 1 SD)
uor_glm_perSD <- function(outcome, zvec, data_mask = NULL) {
  df <- data.frame(y = outcome, z = zvec)
  if (!is.null(data_mask)) df <- df[data_mask, , drop = FALSE]
  fit <- glm(y ~ z, family = binomial, data = df)
  cd  <- summary(fit)$coefficients
  beta <- cd["z", "Estimate"]; se <- cd["z", "Std. Error"]
  list(ORsd = exp(beta), Lsd = exp(beta - 1.96*se), Usd = exp(beta + 1.96*se))
}

## Back-translation helper: raise per-SD OR triplet to (unit / SD_raw)
bt_unit <- function(or_sd_triplet, raw_sd, unit_scale = 1) {
  pow <- unit_scale / raw_sd
  or_sd_triplet^pow
}

## Gated SDs (recreate your gating domains)
sd_gplus  <- function(POut) {
  Gplus <- pmax(POut$Gravida - 1, 0)
  sd(Gplus[POut$Gravida >= 2], na.rm = TRUE)
}
sd_hxhome <- function(POut) {
  sd(POut$hx_home_sum[POut$Gravida >= 2], na.rm = TRUE)
}

## Pull AOR (per SD or vs ref) from GLMM coefficient row
get_aor_ci <- function(glmm, term_regex_exact) {
  cd <- coef_df(glmm)
  j  <- which(grepl(paste0("^", term_regex_exact, "$"), rownames(cd)))
  if (!length(j)) return(NA_character_)
  fmt_or(cd$Estimate[j], cd$`Std. Error`[j])
}

## Pull AOR per-SD triplet (numeric) from GLMM
get_aor_triplet <- function(glmm, term_regex_exact) {
  cd <- coef_df(glmm)
  j  <- which(grepl(paste0("^", term_regex_exact, "$"), rownames(cd)))
  if (!length(j)) return(c(NA,NA,NA))
  beta <- cd$Estimate[j]; se <- cd$`Std. Error`[j]
  c(exp(beta), exp(beta - 1.96*se), exp(beta + 1.96*se))
}

## ---------- Ensure expected factor codings ----------
if (!is.factor(dat2$home_birthYN)) {
  dat2$home_birthYN <- factor(dat2$home_birthYN, levels = c(0,1), labels = c("No","Yes"))
}
## Ensure G1 is a factor for printing unless you already re-coded:
if (!is.null(dat2$G1) && !is.factor(dat2$G1)) {
  dat2$G1 <- factor(dat2$G1, levels = c(0,1), labels = c("Multigravida","Primigravida"))
}

## ---------- Which variables to include ----------
## Categorical (reference = first factor level as encoded in dat2)
cats <- list(
  list(var = "hi2",        label = "Insurance",
       level_map = c("no_insur" = "no insurance", "insur" = "have insurance")),
  list(var = "Age_bin",    label = "Age group",     level_map = NULL),
  list(var = "MEthnicCAT", label = "Ethnicity",     level_map = NULL),
  list(var = "yrproject",  label = "Project year",  level_map = NULL),
  list(var = "G1",         label = "Primigravida indicator",
       level_map = c("Multigravida" = "Multigravida", "Primigravida" = "Primigravida")),
  list(var = "CanRead",    label = "Literacy",      level_map = c("1" = "Literate", "2" = "Not literate"))
)

## Continuous (term = name in GLMM; raw_var = name in POut or "GATED")
## unit_scale is the increment you want the per-unit back-translation for
conts <- list(
  list(term = "NOC_MX_max", label = "Number of ANC visits",
       raw_var = "NOC_MX_max", unit_label = "per visit", unit_scale = 1),
  list(term = "EGA_birth",  label = "Gestational age at birth",
       raw_var = "EGA_birth",  unit_label = "per week",  unit_scale = 1),
  list(term = "HowLongYear",label = "Years in area",
       raw_var = "HowLongYear",unit_label = "per year",  unit_scale = 1),
  list(term = "EnterDays",  label = "Calendar time (in days)",
       raw_var = "EnterDays",  unit_label = "per 30 days", unit_scale = 30),
  list(term = "Gplus_z",    label = "Additional pregnancies (multigravidae)",
       raw_var = "GATED",      unit_label = "per pregnancy", unit_scale = 1),
  list(term = "hxhome_M_z", label = "Prior home-birth history (multigravidae)",
       raw_var = "GATED",      unit_label = "per prior homebirth", unit_scale = 1)
)

## ---------- Build the table ----------
rows <- list()

## A) CATEGORICAL: counts + crude UOR (2×2) + AOR (GLMM)
for (spec in cats) {
  var <- spec$var; lab <- spec$label; m <- spec$level_map
  if (!var %in% names(dat2)) next

  if (!is.factor(dat2[[var]])) dat2[[var]] <- factor(dat2[[var]])
  levs <- levels(droplevels(dat2[[var]]))

  ## Counts
  cnt  <- level_counts(dat2, var, outcome_var = "home_birthYN", yes_label = "Yes")
  mm   <- match(levs, cnt$Level)
  cnt  <- cnt[mm, , drop = FALSE]; row.names(cnt) <- NULL

  ## Reference level counts
  ref_yes <- cnt$Yes[1]; ref_tot <- cnt$Total[1]

  ## Loop levels in order (ref first)
  for (i in seq_along(levs)) {
    L <- levs[i]
    disp_level <- if (!is.null(m) && L %in% names(m)) m[[L]] else L

    if (i == 1L) {
      crude <- "Comparator"
      aor   <- "Comparator"
    } else {
      ## crude UOR from 2×2 counts vs reference
      crude <- or_2x2(yL = cnt$Yes[i], nL = cnt$Total[i], yR = ref_yes, nR = ref_tot)
      ## AOR from GLMM coefficient term "varLevel"
      term_regex <- paste0(esc_re(var), esc_re(L))
      aor <- get_aor_ci(logit_mod_gated, term_regex)
      if (is.na(aor)) aor <- "—"
    }

    rows[[length(rows) + 1]] <- data.frame(
      Variable = lab,
      Category = disp_level,
      Yes      = cnt$Yes[i],
      Total    = cnt$Total[i],
      Percent  = ifelse(is.na(cnt$Percent[i]), NA, round(cnt$Percent[i], 1)),
      `UOR per SD / vs ref` = crude,
      `AOR per SD / vs ref` = aor,
      `UOR per unit`        = "",
      `AOR per unit`        = "",
      stringsAsFactors = FALSE, check.names = FALSE
    )
  }
}

## B) CONTINUOUS: crude per SD (glm), AOR per SD (GLMM), back-translate to per-unit (both)
for (spec in conts) {
  term <- spec$term; lab <- spec$label; raw_name <- spec$raw_var
  unit_lab <- spec$unit_label; unit_scale <- spec$unit_scale

  ## GLMM per SD (from your fitted model)
  a_trip <- get_aor_triplet(logit_mod_gated, esc_re(term))
  a_fmt  <- if (all(is.finite(a_trip))) sprintf("%.2f (%.2f–%.2f)", a_trip[1], a_trip[2], a_trip[3]) else "—"





  ## Raw SD for back-translation & crude per-SD UOR (glm)
  if (raw_name == "GATED") {
    ## ----- Multigravidae-only crude UORs (subset & standardize within MG) -----
    mg <- dat2$Gravida_raw >= 2 & !is.na(dat2$home_birthYN)

    if (term == "Gplus_z") {
      raw_vec <- pmax(POut$Gravida - 1, 0)                 # pregnancies beyond first
      s_raw   <- sd_gplus(POut)                            # SD among multigravidae (for back-translation)
    } else if (term == "hxhome_M_z") {
      raw_vec <- POut$hx_home_sum                          # prior home births
      s_raw   <- sd_hxhome(POut)                           # SD among multigravidae (for back-translation)
    } else {
      next
    }

    ## Standardize within multigravidae only (crude UOR per SD)
    z_mg <- as.numeric(scale(raw_vec[mg]))
    u_sd <- tryCatch(uor_glm_perSD(dat2$home_birthYN[mg], z_mg), error = function(e) NULL)

    if (!is.null(u_sd)) {
      u_fmt  <- sprintf("%.2f (%.2f–%.2f)", u_sd$ORsd, u_sd$Lsd, u_sd$Usd)
      u_trip <- c(u_sd$ORsd, u_sd$Lsd, u_sd$Usd)
    } else {
      u_fmt  <- "—"
      u_trip <- c(NA, NA, NA)
    }

    ## Also compute the GLMM per-SD AOR triplet (already done above as a_trip)
    ## s_raw already set for back-translation

    ## OPTIONAL: annotate label to make subset explicit
    lab_out <- paste0(lab, " (multigravidae only)")

  } else {
    ## ----- Non-gated continuous: crude UOR per SD on full sample -----
    raw   <- POut[[raw_name]]
    s_raw <- sd(raw, na.rm = TRUE)                         # raw SD for back-translation
    z     <- as.numeric(scale(raw))                        # standardize on full sample

    u_sd <- tryCatch(uor_glm_perSD(dat2$home_birthYN, z), error = function(e) NULL)
    if (!is.null(u_sd)) {
      u_fmt  <- sprintf("%.2f (%.2f–%.2f)", u_sd$ORsd, u_sd$Lsd, u_sd$Usd)
      u_trip <- c(u_sd$ORsd, u_sd$Lsd, u_sd$Usd)
    } else {
      u_fmt  <- "—"
      u_trip <- c(NA, NA, NA)
    }

    lab_out <- lab
  }




  ## Back-translate per-unit (both crude and adjusted), if SD is valid
  if (is.finite(s_raw) && s_raw > 0) {
    ## crude per unit
    if (all(is.finite(u_trip))) {
      u_per <- bt_unit(u_trip, s_raw, unit_scale)
      u_unit_fmt <- sprintf("%.2f (%.2f–%.2f)", u_per[1], u_per[2], u_per[3])
    } else {
      u_unit_fmt <- "—"
    }
    ## adjusted per unit
    if (all(is.finite(a_trip))) {
      a_per <- bt_unit(a_trip, s_raw, unit_scale)
      a_unit_fmt <- sprintf("%.2f (%.2f–%.2f)", a_per[1], a_per[2], a_per[3])
    } else {
      a_unit_fmt <- "—"
    }
  } else {
    u_unit_fmt <- "—"; a_unit_fmt <- "—"
  }

  rows[[length(rows) + 1]] <- data.frame(
    Variable = lab_out,
    Category = paste0("per 1 SD  |  per-unit (", unit_lab, ")"),
    Yes = NA, Total = NA, Percent = NA,
    `UOR per SD / vs ref` = u_fmt,
    `AOR per SD / vs ref` = a_fmt,
    `UOR per unit`        = u_unit_fmt,
    `AOR per unit`        = a_unit_fmt,
    stringsAsFactors = FALSE, check.names = FALSE
  )
}

## ---------- Order similar to your narrative ----------
ord_key <- function(v) {
  if (v == "Insurance") return(1)
  if (v == "Age group") return(2)
  if (v == "Ethnicity") return(3)
  if (v == "Project year") return(4)
  if (v == "Primigravida indicator") return(5)
  if (v == "Literacy") return(6)
  if (grepl("^Additional pregnancies", v)) return(7)
  if (grepl("^Prior home-birth history", v)) return(8)
  if (v == "Years in area") return(9)
  if (v == "Calendar time (in days)") return(10)
  if (v == "Gestational age at birth") return(11)
  if (v == "Number of ANC visits") return(12)
  99
}

res <- do.call(rbind, rows)
res$.__k <- vapply(res$Variable, ord_key, numeric(1))
res <- res[order(res$.__k, res$Variable, res$Category), ]
res$.__k <- NULL
row.names(res) <- NULL

## ---------- Output ----------
print(res, row.names = FALSE)
## Optionally save:
# write.csv(res, "overall_results_crude_vs_adjusted_with_perunit.csv", row.names = FALSE)






###################################
#### table for stratified models (25+ age bins; anc4 → NOC_MX_max continuous)

## ===== Side-by-side AOR table WITH reference rows shown =====
## - Shows one row per factor level (incl. the reference: 1.00 (ref))
## - Continuous terms shown once (per SD), no comparator row
## - Works with: logit_mod_prim_3, logit_mod_multi_3
##   (and will adapt if you used different object names)

## ---------- Pick model objects ----------
prim_fit  <- if (exists("logit_mod_prim_3"))  logit_mod_prim_3  else logit_mod_prim
multi_fit <- if (exists("logit_mod_multi_3")) logit_mod_multi_3 else logit_mod_multi

## You must pass the corresponding data frames used to fit:
prim_data  <- dat_prim
multi_data <- dat_multi

## ---------- Helpers ----------
fmt_or <- function(est, se) {
  OR  <- exp(est); LCL <- exp(est - 1.96 * se); UCL <- exp(est + 1.96 * se)
  sprintf("%.2f (%.2f–%.2f)", OR, LCL, UCL)
}
fmt_ref <- function() "1.00 (ref)"
fmt_na  <- function() "—"

esc_re <- function(x) gsub("([][{}()+*^$.|?\\\\])", "\\\\\\1", x)

coef_df <- function(fit) {
  as.data.frame(summary(fit)$coefficients, stringsAsFactors = FALSE)
}

## Find term row for a factor 'var' at level 'lev'
find_term_row <- function(cd, var, lev) {
  pat <- paste0("^", esc_re(var), esc_re(lev), "$")
  which(grepl(pat, rownames(cd)))
}

## Build rows for a factor with explicit reference
build_factor_rows <- function(fit, data, var, prefix_label, level_map = NULL) {
  if (!var %in% names(data)) return(NULL)
  levs <- levels(droplevels(data[[var]]))
  if (is.null(levs)) return(NULL)
  cd <- coef_df(fit)

  ## Map raw level names to pretty labels if provided
  lab <- function(x) {
    if (!is.null(level_map) && x %in% names(level_map)) return(level_map[[x]])
    x
  }

  rows <- lapply(seq_along(levs), function(i) {
    L <- levs[i]
    lbl <- paste0(prefix_label, ": ", lab(L))
    if (i == 1L) {
      data.frame(Label = paste0(lbl, " (ref)"),
                 AOR95 = fmt_ref(), stringsAsFactors = FALSE)
    } else {
      j <- find_term_row(cd, var, L)
      if (length(j) == 1L) {
        data.frame(Label = lbl,
                   AOR95 = fmt_or(cd$Estimate[j], cd$`Std. Error`[j]),
                   stringsAsFactors = FALSE)
      } else {
        data.frame(Label = lbl, AOR95 = fmt_na(), stringsAsFactors = FALSE)
      }
    }
  })
  do.call(rbind, rows)
}

## Binary 0/1 row-pair (e.g., Caesarean01)
build_binary01_rows <- function(fit, var, prefix_label, yes_label = "Yes", no_label = "No") {
  cd <- coef_df(fit)
  if (!any(grepl(paste0("^", esc_re(var), "$"), rownames(cd)))) return(NULL)
  j <- which(grepl(paste0("^", esc_re(var), "$"), rownames(cd)))
  yes <- fmt_or(cd$Estimate[j], cd$`Std. Error`[j])
  data.frame(
    Label = c(paste0(prefix_label, ": ", no_label, " (ref)"),
              paste0(prefix_label, ": ", yes_label)),
    AOR95 = c(fmt_ref(), yes),
    stringsAsFactors = FALSE
  )
}

## Continuous single-row
build_cont_row <- function(fit, var, label) {
  cd <- coef_df(fit)
  j <- which(grepl(paste0("^", esc_re(var), "$"), rownames(cd)))
  if (length(j) != 1L) return(NULL)
  data.frame(Label = label, AOR95 = fmt_or(cd$Estimate[j], cd$`Std. Error`[j]),
             stringsAsFactors = FALSE)
}

## Pretty level mappings you can tweak
map_hi2  <- c("no_insur" = "None", "insur" = "Insured")
## If your CanRead levels are "1"=Yes, "2"=No, use:
map_read <- c("1" = "Yes", "2" = "No")
## If MEthnicCAT has readable labels already, skip mapping.

## ---------- Build model-specific tables ----------
build_model_table <- function(fit, data, which_model = c("prim","multi")) {
  which_model <- match.arg(which_model)

  out <- list()

  ## Factors (all levels incl. ref)
  if ("hi2"        %in% names(data)) out$hi2  <- build_factor_rows(fit, data, "hi2",        "Insurance",    map_hi2)

  ## Age for these stratified models uses the 25+ top bin
  ## Expect variable name: Age_bin25p with levels c("25+","15–19","20–24") (25+ is ref)
  if ("Age_bin25p" %in% names(data)) out$age <- build_factor_rows(fit, data, "Age_bin25p", "Age")

  if ("yrproject"  %in% names(data)) out$yr   <- build_factor_rows(fit, data, "yrproject",  "Project year")
  if ("MEthnicCAT" %in% names(data)) out$eth  <- build_factor_rows(fit, data, "MEthnicCAT", "Ethnicity")
  if ("CanRead"    %in% names(data)) out$read <- build_factor_rows(fit, data, "CanRead",    "Can read", map_read)

  ## Gravidity-related & other continuous
  if (which_model == "multi") {
    out$gplus  <- build_cont_row(fit, "Gplus_z",    "Pregnancies beyond first (per SD)")
    out$csec   <- build_binary01_rows(fit, "Caesarean01", "Prior Caesarean")
    out$hxhome <- build_cont_row(fit, "hxhome_M_z", "Prior home births (per SD)")
  }

  out$howlong <- build_cont_row(fit, "HowLongYear", "Years in area (per SD)")
  out$time    <- build_cont_row(fit, "EnterDays",   "Calendar time (per SD)")
  out$ega     <- build_cont_row(fit, "EGA_birth",   "Gestational age at birth (per SD)")

  ## NEW: replace anc4 (categorical) with NOC_MX_max (continuous)
  if ("NOC_MX_max" %in% names(data)) out$noc <- build_cont_row(fit, "NOC_MX_max", "NOC_MX_max (per SD)")

  do.call(rbind, out)
}

prim_tbl  <- build_model_table(prim_fit,  prim_data,  "prim")
multi_tbl <- build_model_table(multi_fit, multi_data, "multi")

## ---------- Merge side-by-side ----------
all_labels <- unique(c(prim_tbl$Label, multi_tbl$Label))
res <- data.frame(
  Label = all_labels,
  check.names = FALSE, stringsAsFactors = FALSE
)

m1 <- match(res$Label, prim_tbl$Label)
m2 <- match(res$Label, multi_tbl$Label)

res[,"Primigravid AOR (95% CI)"] <- ifelse(is.na(m1), "—", prim_tbl$AOR95[m1])
res[,"Multigravid AOR (95% CI)"] <- ifelse(is.na(m2), "—", multi_tbl$AOR95[m2])

## ---------- Order rows nicely ----------
order_key <- function(lbl) {
  if (grepl("^Insurance:", lbl))               return(1)
  if (grepl("^Age:", lbl))                     return(2)
  if (grepl("^Project year", lbl))             return(3)
  if (grepl("^Ethnicity:", lbl))               return(4)
  if (grepl("^Can read", lbl))                 return(5)
  if (grepl("^Pregnancies beyond first", lbl)) return(6)
  if (grepl("^Prior Caesarean", lbl))          return(7)
  if (grepl("^Prior home births", lbl))        return(8)
  if (grepl("^Years in area", lbl))            return(9)
  if (grepl("^Calendar time", lbl))            return(10)
  if (grepl("^Gestational age", lbl))          return(11)
  if (grepl("^NOC_MX_max", lbl))               return(12)
  99
}

## Within-age block ordering (25+ is reference)
age_ord <- function(lbl) {
  if (!grepl("^Age:", lbl)) return(99)
  if (grepl("25\\+.*\\(ref\\)$", lbl)) return(0)  # ensure ref first
  if (grepl("15–19", lbl))             return(1)
  if (grepl("20–24", lbl))             return(2)
  if (grepl("25\\+", lbl))             return(3)  # non-ref label, if ever present
  50
}

## Keep reference rows first within their blocks
is_ref <- function(lbl) grepl("\\(ref\\)$", lbl)

ord  <- vapply(res$Label, order_key, numeric(1))
aord <- vapply(res$Label, age_ord,   numeric(1))
rord <- ifelse(is_ref(res$Label), -1, 0)

res <- res[order(ord, aord, rord, res$Label), ]
row.names(res) <- NULL

## ---------- Print / save ----------
print(res, row.names = FALSE)
# write.csv(res, "stratified_AOR_table_with_refs_25plus.csv", row.names = FALSE)


write.csv(res, "C:/Users/salem/Dropbox/Rose_Outreach/stratified_AORs.csv", row.names = FALSE)
C:/Users/salem/Dropbox/Rose_Outreach/


out <- "C:/Users/salem/Dropbox/Rose_Outreach/stratified_AORs2.xlsx"
writexl::write_xlsx(
  list(
    "Table1_Descripts_UOR_AOR"     = res
  ),
  path = out
)




##############################################################################
########## probability curves for ANC visits  — with marginalization over u
##############################################################################
## fixed-effects draws + integration over the random intercept distribution
## optional sigma draws via profile CI (fallback: hold sigma fixed)
##############################################################################

library(ggplot2)
library(scales)
library(mvtnorm)
library(lme4)

## ---- Setup for prediction on RAW ANC while model uses Z-scores ----

## 1) Get the exact rows used in the model and their raw ANC values
mf_fit <- model.frame(logit_mod_gated)              # data actually used by the model
row_idx <- as.integer(rownames(mf_fit))             # should map to POut rows
anc_raw_used <- POut$NOC_MX_max[row_idx]

## 2) Centering/scaling parameters to go from RAW -> Z exactly as in dat2
anc_mean <- mean(anc_raw_used, na.rm = TRUE)
anc_sd   <- sd(anc_raw_used,   na.rm = TRUE)

## 3) Build the RAW ANC grid from observed *raw* range (integers only)
anc_grid_raw <- seq(
  floor(min(anc_raw_used, na.rm = TRUE)),
  ceiling(max(anc_raw_used, na.rm = TRUE)),
  by = 1
)
## Convert RAW -> Z for the model
anc_grid_z <- (anc_grid_raw - anc_mean) / anc_sd

## 4) Base frames for marginal standardization (hold cohort mix as observed)
base_all   <- mf_fit
base_primi <- mf_fit[mf_fit$G1 == 1, , drop = FALSE]
base_multi <- mf_fit[mf_fit$G1 == 0, , drop = FALSE]

## Helper to build model matrix that matches fitted design columns
build_X <- function(model, newdata){
  mm <- model.matrix(terms(model), newdata)
  Xc <- colnames(getME(model, "X"))
  miss <- setdiff(Xc, colnames(mm))
  if (length(miss))
    mm <- cbind(mm, matrix(0, nrow(mm), length(miss), dimnames = list(NULL, miss)))
  mm[, Xc, drop = FALSE]
}
make_X_list <- function(base_df, zgrid, model){
  lapply(zgrid, function(z){
    nd <- base_df
    nd$NOC_MX_max <- z             # IMPORTANT: model expects z-scores internally
    build_X(model, nd)
  })
}

X_overall <- make_X_list(base_all,   anc_grid_z, logit_mod_gated)
X_primi   <- make_X_list(base_primi, anc_grid_z, logit_mod_gated)
X_multi   <- make_X_list(base_multi, anc_grid_z, logit_mod_gated)

## ---- Marginalization over random intercept u ~ N(0, sigma^2) ----
## Compute E_u[logit^{-1}(eta + u)] by Monte Carlo; reuse u draws for speed
marginalize_over_u <- function(eta_vec, sd_u, K = 200L){
  if (sd_u <= 0) return(mean(plogis(eta_vec)))
  u <- rnorm(K, 0, sd_u)
  ## broadcast eta across K draws, then average over u and over individuals
  p <- plogis(matrix(eta_vec, nrow = length(eta_vec), ncol = K) +
                rep(u, each = length(eta_vec)))
  mean(p)
}

## Fast point estimates (population-averaged): integrate over u at the MLEs
sigma_hat <- as.numeric(attr(VarCorr(logit_mod_gated)[[1]], "stddev"))
beta_hat  <- fixef(logit_mod_gated)

point_curve_from_eta <- function(X_list, beta, sd_u, K = 400L){
  vapply(X_list, function(Xk){
    eta <- as.vector(Xk %*% beta)
    marginalize_over_u(eta, sd_u, K = K)
  }, numeric(1))
}

fig1_overall <- point_curve_from_eta(X_overall, beta_hat, sigma_hat)
fig1_primi   <- point_curve_from_eta(X_primi,   beta_hat, sigma_hat)
fig1_multi   <- point_curve_from_eta(X_multi,   beta_hat, sigma_hat)

## ---- Uncertainty: draws of beta and (optionally) sigma ----
## Draw betas from MVN; try to obtain sigma draws from profile CI on SD
set.seed(123)
B <- 600L
Vb <- as.matrix(vcov(logit_mod_gated))
beta_sims <- rmvnorm(B, beta_hat, Vb)

## Try to get profile CI for the random-intercept SD (theta on SD scale in lme4)
## If profiling is unavailable/slow, fall back to fixed sigma
get_sigma_draws <- function(mod, B, sigma_hat){
  sigma_sims <- rep(sigma_hat, B)
  ok <- FALSE
  ## Attempt quick profile; if unavailable, skip
  suppressWarnings({
    p <- try(confint(mod, parm = "theta_", method = "profile", oldNames = FALSE), silent = TRUE)
  })
  if (!inherits(p, "try-error")) {
    ## row name typically "sd_(Intercept)|<group>"
    srow <- grep("^sd_\\(Intercept\\)\\|", rownames(p), value = TRUE)
    if (length(srow) == 1L) {
      lo <- p[srow, 1]; hi <- p[srow, 2]
      if (is.finite(lo) && lo > 0 && is.finite(hi) && hi > 0 && hi > lo) {
        ## Approximate sampling for SD via log-normal fitted to CI endpoints
        sdlog <- (log(hi) - log(lo)) / (2 * 1.96)
        meanlog <- log(sigma_hat)
        sigma_sims <- rlnorm(B, meanlog = meanlog, sdlog = sdlog)
        ok <- TRUE
      }
    }
  }
  attr(sigma_sims, "profile_based") <- ok
  sigma_sims
}

sigma_sims <- get_sigma_draws(logit_mod_gated, B, sigma_hat)

## For each draw, integrate over u and average across individuals (marginal standardization)
mean_probs_from_draws <- function(X_list, betas, sigma_draws, K_u = 200L){
  B <- nrow(betas); K <- length(X_list)
  out <- matrix(NA_real_, B, K)
  for (k in seq_len(K)){
    Xk <- X_list[[k]]
    for (b in seq_len(B)){
      eta <- as.vector(Xk %*% betas[b, ])
      out[b, k] <- marginalize_over_u(eta, sigma_draws[b], K = K_u)
    }
  }
  out
}

sim_overall <- mean_probs_from_draws(X_overall, beta_sims, sigma_sims, K_u = 200L)
sim_primi   <- mean_probs_from_draws(X_primi,   beta_sims, sigma_sims, K_u = 200L)
sim_multi   <- mean_probs_from_draws(X_multi,   beta_sims, sigma_sims, K_u = 200L)

ci2 <- function(M)
  rbind(lwr = apply(M, 2, quantile, 0.025),
        upr = apply(M, 2, quantile, 0.975))

## 7) Assemble Figure 1 data (RAW x-axis, z used internally for preds)
pred_fig1 <- rbind(
  data.frame(ANC = anc_grid_raw, Group = "Overall",
             Pred = fig1_overall, t(ci2(sim_overall))),
  data.frame(ANC = anc_grid_raw, Group = "Primigravida",
             Pred = fig1_primi,   t(ci2(sim_primi))),
  data.frame(ANC = anc_grid_raw, Group = "Multigravida",
             Pred = fig1_multi,   t(ci2(sim_multi)))
)

## 8) Figure 2: Multigravidas, prior C-section vs none (toggle only Csec_M)
make_X_list_csec <- function(base_df, zgrid, model, csec_val){
  lapply(zgrid, function(z){
    nd <- base_df
    nd$NOC_MX_max <- z
    nd$Csec_M     <- csec_val
    build_X(model, nd)
  })
}
X_csec0 <- make_X_list_csec(base_multi, anc_grid_z, logit_mod_gated, 0)
X_csec1 <- make_X_list_csec(base_multi, anc_grid_z, logit_mod_gated, 1)

fig2_csec0 <- point_curve_from_eta(X_csec0, beta_hat, sigma_hat)
fig2_csec1 <- point_curve_from_eta(X_csec1, beta_hat, sigma_hat)

sim_csec0 <- mean_probs_from_draws(X_csec0, beta_sims, sigma_sims, K_u = 200L)
sim_csec1 <- mean_probs_from_draws(X_csec1, beta_sims, sigma_sims, K_u = 200L)

pred_fig2 <- rbind(
  data.frame(ANC = anc_grid_raw, Group = "No prior C-section",
             Pred = fig2_csec0, t(ci2(sim_csec0))),
  data.frame(ANC = anc_grid_raw, Group = "Prior C-section",
             Pred = fig2_csec1, t(ci2(sim_csec1)))
)

## 9) Plot helpers (integer ticks over observed RAW ANC)
plot_curve <- function(df, title, subtitle){
  ggplot(df, aes(x = ANC, y = Pred, color = Group, fill = Group)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.18, color = NA) +
    geom_line(size = 0.7) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0,1)) +
    scale_x_continuous(breaks = sort(unique(df$ANC))) +
    labs(x = "Number of ANC visits",
         y = "Predicted probability of home birth",
         title = title, subtitle = subtitle) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top")
}

fig1 <- plot_curve(
  pred_fig1,
  "",
  ""
  # "Home-birth probability by ANC visits (marginal standardized)",
  # "Population-averaged predictions (GLMM): uncertainty from fixed effects and random-intercept variance; integrates over u ~ N(0, σ²)."
)
fig2 <- plot_curve(
  pred_fig2,
  "",
  ""
  # "Home-birth probability among multigravidas by prior C-section",
  # "Population-averaged predictions (GLMM): uncertainty from fixed effects and random-intercept variance; integrates over u ~ N(0, σ²)."
)

print(fig1); print(fig2)

## Optional: message about sigma draws source
if (isTRUE(attr(sigma_sims, "profile_based"))) {
  message("CIs reflect uncertainty in beta and sigma (sigma via profile-based log-normal draws).")
} else {
  message("CIs reflect uncertainty in beta; sigma held at MLE (no profile CI available).")
}
## Optional: save publication PNGs (adjust width/height as needed)
#ggsave("C:/Users/salem/Dropbox/Rose_Outreach/fig1_overall_by_grav.png", fig1, width = 6.5, height = 4.2, dpi = 300)
ggsave("C:/Users/salem/Dropbox/Rose_Outreach/fig2_multigrav_csection.png", fig2, width = 6.5, height = 4.2, dpi = 300)


## ---- Figure 1b: Primi vs Multigravida only (no Overall) ----

## Primigravida vs Multigravida only, using existing objects/functions:
pred_fig1_primi_multi <- subset(pred_fig1, Group %in% c("Primigravida","Multigravida"))

## (Optional) ensure ordering in legend:
pred_fig1_primi_multi$Group <- factor(
  pred_fig1_primi_multi$Group,
  levels = c("Primigravida","Multigravida")
)

fig1b <- plot_curve(
  pred_fig1_primi_multi,
  "",
  ""
  # "Home-birth probability by ANC visits (marginal standardized)",
  # "Population-averaged predictions (GLMM): uncertainty from fixed effects and random-intercept variance; integrates over u ~ N(0, σ²)."
)

print(fig1b)

## Optional: save to file
ggsave("C:/Users/salem/Dropbox/Rose_Outreach/fig1b_primi_vs_multi.png", fig1b, width = 6.5, height = 4.2, dpi = 300)









## 10) Quick sanity checks you asked for
cat("\nObserved RAW ANC range (model rows):\n"); print(range(anc_raw_used, na.rm = TRUE))
cat("\nCrude Yes/Total at ANC = 1, 4, 8 (overall):\n")
for (k in c(1,4,8)) {
  ok <- POut$NOC_MX_max[row_idx] == k
  yes <- sum(dat2$home_birthYN[row_idx][ok] == "Yes", na.rm = TRUE)
  tot <- sum(ok, na.rm = TRUE)
  cat(sprintf("ANC=%d: %d/%d (%.1f%%)\n", k, yes, tot, 100*yes/pmax(tot,1)))
}
























