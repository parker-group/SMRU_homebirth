
# -----------------------
# Load packages
# -----------------------


library(lme4)
library(boot)
library(ggplot2)
library(scales)
library(mvtnorm)

# -------------------------------
# 0) Define key paths + ensure output folders exist
# -------------------------------
DATA_FILE <- "data/homebirth_analysis.csv"
dir.create("output", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# -------------------------------
# 1) Load dataset
# -------------------------------
POut <- read.csv(DATA_FILE, header = TRUE, sep = ",")

# Rename pseudonymized IDs to the canonical names used throughout the script
names(POut)[names(POut) == "ANC_id"]     <- "ANC_code"
names(POut)[names(POut) == "OR_site_id"] <- "OR_site"

# Safety checks: required columns present + no missing IDs + unique person rows
if (!("ANC_code" %in% names(POut))) stop("Missing ANC_code in data/homebirth_analysis.csv")
if (!("OR_site"  %in% names(POut))) stop("Missing OR_site in data/homebirth_analysis.csv")

if (anyNA(POut$ANC_code)) stop("ANC_code has missing values in data/homebirth_analysis.csv")
if (anyNA(POut$OR_site))  stop("OR_site has missing values in data/homebirth_analysis.csv")

if (anyDuplicated(POut$ANC_code)) {
  dup <- unique(POut$ANC_code[duplicated(POut$ANC_code)])
  stop("Expected 1 row per person, but ANC_code is duplicated (examples): ",
       paste(head(dup), collapse = ", "))
}

# Treat IDs as IDs (prevents accidental numeric behavior)
POut$ANC_code <- factor(POut$ANC_code)
POut$OR_site  <- factor(POut$OR_site)

# -------------------------------
# 2) Ensure that dates are seen as dates, and generate a day on which a participant entered the study
# -------------------------------
POut$EnterDate <- as.Date(POut$EnterDate, format = "%m/%d/%Y")
POut$EnterDays <- as.numeric(POut$EnterDate - min(POut$EnterDate, na.rm = TRUE))

# -------------------------------
# 3) Build analysis frame from the public dataset
# -------------------------------
vars_keep <- c(
  "ANC_code","OR_site",
  "home_birthYN",
  "hi","Age","yrproject","MEthnicCAT","CanRead",
  "Parity","Caesarean","hx_home_sum",
  "HowLongYear","EnterDate","EnterDays",
  "enrolEGA","EGA_birth",
  "NOC_MX_max",
  "prematureY_N"
)

missing <- setdiff(vars_keep, names(POut))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse = ", "))

dat2 <- POut[, vars_keep]

# -------------------------------
# 4) Factor coding and age bin system (Age_bin25p)
# -------------------------------
dat2$home_birthYN <- factor(dat2$home_birthYN, levels = c(0,1), labels = c("No","Yes"))

dat2$MEthnicCAT <- factor(dat2$MEthnicCAT)
dat2$OR_site    <- factor(dat2$OR_site)
dat2$yrproject  <- factor(dat2$yrproject)

# Keep hi (health insurance) as a factor; ensure it's character-safe first
dat2$hi <- as.character(dat2$hi)
dat2$hi <- factor(dat2$hi)

dat2$Caesarean <- factor(dat2$Caesarean)
dat2$CanRead   <- factor(dat2$CanRead, levels = c(0,1), labels = c("No","Yes"))
dat2$prematureY_N <- factor(dat2$prematureY_N)

# Insurance collapse (hi2): insur vs no_insur (ref = no_insur)
dat2$hi2 <- ifelse(dat2$hi %in% c("1","2"), "insur", "no_insur")
dat2$hi2 <- factor(dat2$hi2, levels = c("no_insur","insur"))

# SINGLE age system everywhere: 15–19, 20–24, 25+ ; ref = 25+
dat2$Age_bin25p <- cut(dat2$Age,
                       breaks = c(15, 20, 25, Inf),
                       right  = FALSE,
                       labels = c("15–19","20–24","25+"))
dat2$Age_bin25p <- factor(dat2$Age_bin25p, levels = c("25+","15–19","20–24"))

# -------------------------------
# 4.1) Late first ANC visit indicator (based on RAW enrolEGA weeks)
#     late_first_anc = Yes if enrolEGA > 28 weeks
#     Ref level = No (<=28)
# -------------------------------
dat2$late_first_anc <- ifelse(is.na(dat2$enrolEGA), NA_character_,
                              ifelse(dat2$enrolEGA > 28, "Yes", "No"))
dat2$late_first_anc <- factor(dat2$late_first_anc, levels = c("No","Yes"))

# -------------------------------
# 5) Standardize continuous covariates (mean 0, SD 1) for model stability
#     
# -------------------------------
safe_scale <- function(x) {
  if (all(is.na(x))) return(x)
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(x*0)
  as.numeric((x - m)/s)
}

cont_vars <- c("hx_home_sum","HowLongYear","EnterDays","NOC_MX_max","enrolEGA")
dat2[cont_vars] <- lapply(dat2[cont_vars], safe_scale)

# -------------------------------
# 5.1) Preserve raw counts BEFORE gating (used for multip-only standardized terms)
# -------------------------------
dat2$Parity_raw      <- POut$Parity
dat2$hx_home_sum_raw <- POut$hx_home_sum

# -------------------------------
# 6) Create parity gating variables and multip-only standardized terms
#     - G1: Primip vs Multip (factor)
#     - Gplus_z: births beyond first, z-scored within multiparas
#     - hxhome_M_z: prior home births, z-scored within multiparas
#     - Csec_M: prior C-section, gated to multiparas (factor No/Yes)
# -------------------------------
dat2$G1    <- ifelse(is.na(dat2$Parity_raw), NA_integer_, as.integer(dat2$Parity_raw == 0L))
dat2$M     <- ifelse(is.na(dat2$Parity_raw), NA_integer_, as.integer(dat2$Parity_raw >= 1L))
dat2$Gplus <- ifelse(is.na(dat2$Parity_raw), NA_real_, pmax(dat2$Parity_raw - 1, 0))

dat2$G1 <- factor(ifelse(is.na(dat2$Parity_raw), NA_integer_, as.integer(dat2$Parity_raw == 0L)),
                 levels = c(0, 1), labels = c("Multip","Primip"))

# Robust prior C-section binary (handles inconsistent Yes/No encodings)
dat2$Caesarean01 <- as.integer(as.character(dat2$Caesarean) %in% c("Yes","Ye","1","Y","TRUE","True"))

# Gate prior home births to multiparas; primip gets 0
dat2$hxhome_M <- ifelse(dat2$M == 1L, dat2$hx_home_sum_raw, 0)

# Standardize Gplus and hxhome_M within multiparas only (zeros elsewhere)
dat2$Gplus_z    <- 0
dat2$hxhome_M_z <- 0

idx_M <- which(dat2$M == 1L & !is.na(dat2$Gplus))
if (length(idx_M) > 2) {
  mu_gp <- mean(dat2$Gplus[idx_M], na.rm = TRUE)
  sd_gp <-  sd(dat2$Gplus[idx_M],  na.rm = TRUE)
  if (is.finite(sd_gp) && sd_gp > 0) dat2$Gplus_z[idx_M] <- (dat2$Gplus[idx_M] - mu_gp) / sd_gp
}

idx_H <- which(dat2$M == 1L & !is.na(dat2$hxhome_M))
if (length(idx_H) > 2) {
  mu_hx <- mean(dat2$hxhome_M[idx_H], na.rm = TRUE)
  sd_hx <-  sd(dat2$hxhome_M[idx_H],  na.rm = TRUE)
  if (is.finite(sd_hx) && sd_hx > 0) dat2$hxhome_M_z[idx_H] <- (dat2$hxhome_M[idx_H] - mu_hx) / sd_hx
}

# Prior C-section (gated): multip gets Caesarean01 (NA->0), primip forced 0
dat2$Csec_M <- ifelse(dat2$M == 1L,
                      ifelse(is.na(dat2$Caesarean01), 0L, dat2$Caesarean01),
                      0L)
dat2$Csec_M <- factor(dat2$Csec_M, levels = c(0,1), labels = c("No","Yes"))

# -------------------------------
# 7) Fit pooled GLMM (random intercept by OR_site)
#
# -------------------------------
logit_mod_gated <- glmer(
  home_birthYN ~ hi2 + Age_bin25p + MEthnicCAT + CanRead +
    G1 + Gplus_z + Csec_M + hxhome_M_z + HowLongYear + EnterDays + enrolEGA + NOC_MX_max +
    (1 | OR_site),
  data   = dat2,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# -------------------------------
# 8) OR_site clustering metrics (VAR, MOR, ICC) + bootstrap CIs
#     - Extract random-intercept variance (sigma^2)
#     - Compute MOR and ICC from sigma^2
#     - Bootstrap via bootMer for uncertainty intervals
# -------------------------------

# Function to extract the random-intercept variance (sigma^2) for OR_site
# from the fitted GLMM. This is the between-site variance on the log-odds scale.
get_sigma2 <- function(fit, group = "OR_site") {
  vc <- VarCorr(fit)
  df <- as.data.frame(vc)
  row <- df[df$grp == group & df$var1 == "(Intercept)" & is.na(df$var2), , drop = FALSE]
  if (nrow(row) == 0L) return(NA_real_)
  as.numeric(row$vcov[1])
}


# Function to convert the random-intercept variance into the
# Median Odds Ratio (MOR), which expresses between-site variation
# on the odds-ratio scale for easier interpretation.
mor_from_var <- function(s2) {
  if (!is.finite(s2) || s2 < 0) return(NA_real_)
  exp(0.6744898 * sqrt(2 * s2))
}

# Function to compute the Intraclass Correlation Coefficient (ICC)
# for a logistic mixed model using the standard latent-variable
# approximation where the individual-level residual variance = π² / 3.
icc_from_var <- function(s2) {
  if (!is.finite(s2) || s2 < 0) return(NA_real_)
  s2 / (s2 + (pi^2 / 3))
}

# Helper function to compute percentile-based confidence intervals
# (used later for bootstrap confidence intervals).
pct <- function(x, probs = c(0.025, 0.975)) {
  as.numeric(quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

# Point estimates
sigma2_hat <- get_sigma2(logit_mod_gated, group = "OR_site")
MOR_hat    <- mor_from_var(sigma2_hat)
ICC_hat    <- icc_from_var(sigma2_hat)

# Bootstrap setup (increase B if desired, really should be at least 500)
set.seed(20250908)
B <- 500


# Define the statistic function that will be evaluated during bootstrapping.
# For each bootstrap model refit, we extract:
#   VAR = random-intercept variance
#   MOR = Median Odds Ratio
#   ICC = Intraclass Correlation
stat_fun <- function(fit) {
  s2 <- get_sigma2(fit, group = "OR_site")
  c(VAR = s2,
    MOR = mor_from_var(s2),
    ICC = icc_from_var(s2))
}

# Run the parametric bootstrap using bootMer().
# Each bootstrap draw simulates new outcomes under the fitted model,
# refits the model, and recomputes the clustering statistics.
bs <- bootMer(
  x        = logit_mod_gated,
  FUN      = stat_fun,
  nsim     = B,
  use.u    = FALSE,
  type     = "parametric",
  parallel = "no",
  re.form  = NULL,
  verbose  = 0
)

# Summarize bootstrap uncertainty
SEs    <- apply(bs$t, 2, sd, na.rm = TRUE)
VAR_CI <- pct(bs$t[, "VAR"])
MOR_CI <- pct(bs$t[, "MOR"])
ICC_CI <- pct(bs$t[, "ICC"])

# Assemble all clustering statistics into a single list
# for convenient inspection and reporting.
mor_icc_out <- list(
  sigma2 = unname(sigma2_hat),
  MOR    = unname(MOR_hat),
  ICC    = unname(ICC_hat),
  SE     = setNames(SEs, names(SEs)),
  VAR_CI = setNames(VAR_CI, c("lower","upper")),
  MOR_CI = setNames(MOR_CI, c("lower","upper")),
  ICC_CI = setNames(ICC_CI, c("lower","upper")),
  boot_reps      = B,
  boot_converged = sum(complete.cases(bs$t[, "VAR"]))
)

print(mor_icc_out)

# -------------------------------
# 9) Unified RESULTS TABLE — pooled model
#     Components:
#     (A) Descriptives: Yes/Total/% by category
#     (B) UOR vs ref for categorical predictors (2x2 OR with continuity correction)
#     (C) AOR vs ref for categorical predictors (from pooled GLMM)
#     (D) For continuous predictors:
#         - UOR per 1 SD (simple logistic regression)
#         - AOR per 1 SD (from pooled GLMM)
#         - Back-translate both to per-unit ORs using raw SDs (from POut)
# -------------------------------
coef_df <- function(fit) as.data.frame(summary(fit)$coefficients, stringsAsFactors = FALSE)

fmt_or <- function(beta, se) {
  OR <- exp(beta); L <- exp(beta - 1.96*se); U <- exp(beta + 1.96*se)
  sprintf("%.2f (%.2f–%.2f)", OR, L, U)
}

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

or_2x2 <- function(yL, nL, yR, nR) {
  a <- yL; b <- nL - yL; c <- yR; d <- nR - yR
  if (any(c(a,b,c,d) == 0)) { a <- a + 0.5; b <- b + 0.5; c <- c + 0.5; d <- d + 0.5 }
  OR  <- (a/b)/(c/d)
  seL <- sqrt(1/a + 1/b + 1/c + 1/d)
  L   <- exp(log(OR) - 1.96*seL)
  U   <- exp(log(OR) + 1.96*seL)
  sprintf("%.2f (%.2f–%.2f)", OR, L, U)
}

uor_glm_perSD <- function(outcome, zvec, data_mask = NULL) {
  df <- data.frame(y = outcome, z = zvec)
  if (!is.null(data_mask)) df <- df[data_mask, , drop = FALSE]
  fit <- glm(y ~ z, family = binomial, data = df)
  cd  <- summary(fit)$coefficients
  beta <- cd["z", "Estimate"]; se <- cd["z", "Std. Error"]
  list(ORsd = exp(beta), Lsd = exp(beta - 1.96*se), Usd = exp(beta + 1.96*se))
}

bt_unit <- function(or_sd_triplet, raw_sd, unit_scale = 1) {
  pow <- unit_scale / raw_sd
  or_sd_triplet^pow
}

sd_gplus  <- function(POut) {
  Gplus <- pmax(POut$Parity - 1, 0)
  sd(Gplus[POut$Parity >= 1], na.rm = TRUE)
}
sd_hxhome <- function(POut) {
  sd(POut$hx_home_sum[POut$Parity >= 1], na.rm = TRUE)
}

esc_re  <- function(x) gsub("([][{}()+*^$.|?\\\\])", "\\\\\\1", x)

get_aor_ci <- function(glmm, term_key) {
  cd <- coef_df(glmm)
  j  <- which(rownames(cd) == term_key)
  if (!length(j)) j <- which(grepl(paste0("^", esc_re(term_key), "$"), rownames(cd)))
  if (!length(j)) return(NA_character_)
  fmt_or(cd$Estimate[j], cd$`Std. Error`[j])
}

get_aor_triplet <- function(glmm, term_key) {
  cd <- coef_df(glmm)
  j  <- which(rownames(cd) == term_key)
  if (!length(j)) j <- which(grepl(paste0("^", esc_re(term_key), "$"), rownames(cd)))
  if (!length(j)) return(c(NA, NA, NA))
  beta <- cd$Estimate[j]; se <- cd$`Std. Error`[j]
  c(exp(beta), exp(beta - 1.96*se), exp(beta + 1.96*se))
}

if (!is.factor(dat2$home_birthYN)) {
  dat2$home_birthYN <- factor(dat2$home_birthYN, levels = c("No","Yes"))
}

# (A) Categorical variable specs included in the results table
cats <- list(
  list(var = "hi2",        label = "Insurance",
       level_map = c("no_insur" = "no insurance", "insur" = "have insurance")),
  list(var = "Age_bin25p", label = "Age group",   level_map = NULL),
  list(var = "MEthnicCAT", label = "Ethnicity",   level_map = NULL),
  list(var = "G1",         label = "Primiparous indicator",
       level_map = c("Multip" = "Multip", "Primip" = "Primip")),
  list(var = "CanRead",    label = "Literacy",    level_map = NULL),
  list(var = "Csec_M",     label = "Prior C-section (multiparous gating)",
       level_map = c("No" = "No", "Yes" = "Yes"))
)

# (B) Continuous variable specs: include per-SD and per-unit back-translation
conts <- list(
  list(term = "NOC_MX_max", label = "Number of ANC visits",
       raw_var = "NOC_MX_max", unit_label = "per visit", unit_scale = 1),
  list(term = "enrolEGA",   label = "Gestational age at enrollment",
       raw_var = "enrolEGA", unit_label = "per week",  unit_scale = 1),
  list(term = "HowLongYear",label = "Years in area",
       raw_var = "HowLongYear",unit_label = "per year",  unit_scale = 1),
  list(term = "EnterDays",  label = "Calendar time (in days)",
       raw_var = "EnterDays",  unit_label = "per 30 days", unit_scale = 30),
  list(term = "Gplus_z",    label = "Additional births (multiparous)",
       raw_var = "GATED",      unit_label = "per birth", unit_scale = 1),
  list(term = "hxhome_M_z", label = "Prior home-birth history (multiparous)",
       raw_var = "GATED",      unit_label = "per prior homebirth", unit_scale = 1)
)

rows <- list()

# 9A) Categorical predictors: descriptives + UOR vs ref + AOR vs ref
for (spec in cats) {
  var <- spec$var; lab <- spec$label; m <- spec$level_map
  if (!var %in% names(dat2)) next

  if (!is.factor(dat2[[var]])) dat2[[var]] <- factor(dat2[[var]])
  levs <- levels(droplevels(dat2[[var]]))

  cnt  <- level_counts(dat2, var, outcome_var = "home_birthYN", yes_label = "Yes")
  mm   <- match(levs, cnt$Level)
  cnt  <- cnt[mm, , drop = FALSE]; row.names(cnt) <- NULL

  ref_yes <- cnt$Yes[1]; ref_tot <- cnt$Total[1]

  for (i in seq_along(levs)) {
    L <- levs[i]
    disp_level <- if (!is.null(m) && L %in% names(m)) m[[L]] else L

    if (i == 1L) {
      crude <- "Comparator"
      aor   <- "Comparator"
    } else {
      crude <- or_2x2(yL = cnt$Yes[i], nL = cnt$Total[i], yR = ref_yes, nR = ref_tot)
      term_key <- paste0(var, L)
      aor <- get_aor_ci(logit_mod_gated, term_key)
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

# 9B) Continuous predictors: UOR per SD + AOR per SD + per-unit back-translation
for (spec in conts) {
  term <- spec$term; lab <- spec$label; raw_name <- spec$raw_var
  unit_lab <- spec$unit_label; unit_scale <- spec$unit_scale

  # Adjusted per-SD from pooled model
  a_trip <- get_aor_triplet(logit_mod_gated, term)
  a_fmt  <- if (all(is.finite(a_trip))) sprintf("%.2f (%.2f–%.2f)", a_trip[1], a_trip[2], a_trip[3]) else "—"

  # Unadjusted per-SD from simple logistic regression
  if (raw_name == "GATED") {
    mp <- dat2$Parity_raw >= 1 & !is.na(dat2$home_birthYN)

    if (term == "Gplus_z") {
      raw_vec <- pmax(POut$Parity - 1, 0)
      s_raw   <- sd_gplus(POut)
    } else if (term == "hxhome_M_z") {
      raw_vec <- POut$hx_home_sum
      s_raw   <- sd_hxhome(POut)
    } else next

    z_mp <- as.numeric(scale(raw_vec[mp]))
    u_sd <- tryCatch(uor_glm_perSD(dat2$home_birthYN[mp], z_mp), error = function(e) NULL)

    if (!is.null(u_sd)) {
      u_fmt  <- sprintf("%.2f (%.2f–%.2f)", u_sd$ORsd, u_sd$Lsd, u_sd$Usd)
      u_trip <- c(u_sd$ORsd, u_sd$Lsd, u_sd$Usd)
    } else {
      u_fmt  <- "—"
      u_trip <- c(NA, NA, NA)
    }

    lab_out <- paste0(lab, " (multiparous only)")

  } else {
    # Align crude SD calculation to the exact rows used by the fitted model
    mf_fit  <- model.frame(logit_mod_gated)
    row_idx <- as.integer(rownames(mf_fit))

    y_use   <- dat2$home_birthYN[row_idx]
    raw_use <- POut[[raw_name]][row_idx]

    s_raw <- sd(raw_use, na.rm = TRUE)
    z     <- as.numeric(scale(raw_use))

    u_sd <- tryCatch(uor_glm_perSD(y_use, z), error = function(e) NULL)

    if (!is.null(u_sd)) {
      u_fmt  <- sprintf("%.2f (%.2f–%.2f)", u_sd$ORsd, u_sd$Lsd, u_sd$Usd)
      u_trip <- c(u_sd$ORsd, u_sd$Lsd, u_sd$Usd)
    } else {
      u_fmt  <- "—"
      u_trip <- c(NA, NA, NA)
    }

    lab_out <- lab
  }

  # Back-translate per-SD ORs into per-unit ORs using raw SD in the analysis rows
  if (is.finite(s_raw) && s_raw > 0) {
    if (all(is.finite(u_trip))) {
      u_per <- bt_unit(u_trip, s_raw, unit_scale)
      u_unit_fmt <- sprintf("%.2f (%.2f–%.2f)", u_per[1], u_per[2], u_per[3])
    } else u_unit_fmt <- "—"

    if (all(is.finite(a_trip))) {
      a_per <- bt_unit(a_trip, s_raw, unit_scale)
      a_unit_fmt <- sprintf("%.2f (%.2f–%.2f)", a_per[1], a_per[2], a_per[3])
    } else a_unit_fmt <- "—"
  } else {
    u_unit_fmt <- "—"
    a_unit_fmt <- "—"
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

# 9C) Order and display the pooled results table
ord_key <- function(v) {
  if (v == "Insurance") return(1)
  if (v == "Age group") return(2)
  if (v == "Ethnicity") return(3)
  if (v == "Primiparous indicator") return(4)
  if (v == "Literacy") return(5)
  if (v == "Prior C-section (multiparous gating)") return(6)
  if (grepl("^Additional births", v)) return(7)
  if (grepl("^Prior home-birth history", v)) return(8)
  if (v == "Years in area") return(9)
  if (v == "Calendar time \\(in days\\)") return(10)
  if (v == "Gestational age at enrollment") return(11)
  if (v == "Number of ANC visits") return(12)
  99
}

res <- do.call(rbind, rows)
res$.__k <- vapply(res$Variable, ord_key, numeric(1))
res <- res[order(res$.__k, res$Variable, res$Category), ]
res$.__k <- NULL
row.names(res) <- NULL

print(res, row.names = FALSE)
# write.csv(res, "output/tables/Table1_Descripts_UOR_AOR.csv", row.names = FALSE)

# -------------------------------
# 10) Probability curves for ANC visits (population-averaged over random intercept)
#
#     Approach:
#       - Use the fitted model frame (exact rows used in logit_mod_gated)
#       - Vary raw ANC (x-axis) but set model’s z-scored NOC_MX_max accordingly
#       - Integrate over u ~ N(0, sigma^2) to obtain population-averaged probabilities
#       - Uncertainty via MVN draws of beta (+ optional sigma draws)
# -------------------------------

# 10.1) Pull the exact analysis rows used by the fitted model (ensures perfect alignment)
mf_fit <- model.frame(logit_mod_gated)
# Recover the row indices of those observations in the original dataset.
row_idx <- as.integer(rownames(mf_fit))
# Extract the raw ANC visit counts corresponding to those rows.
anc_raw_used <- POut$NOC_MX_max[row_idx]

# 10.2) Recreate the raw->z mapping for ANC used in the fitted dat2 rows
# The model used a z-scored version of ANC visits, so we compute the
# mean and standard deviation from the model rows.
anc_mean <- mean(anc_raw_used, na.rm = TRUE)
anc_sd   <- sd(anc_raw_used,   na.rm = TRUE)

# Create a grid of ANC visit counts covering the observed range.
# These values form the x-axis of the probability curve.
anc_grid_raw <- seq(
  floor(min(anc_raw_used, na.rm = TRUE)),
  ceiling(max(anc_raw_used, na.rm = TRUE)),
  by = 1
)

# Convert those raw values to the standardized scale used in the model
anc_grid_z <- (anc_grid_raw - anc_mean) / anc_sd

# 10.3) Base frames for marginal standardization (parity subsets)
# Split the model dataset into parity strata.
# These are used to generate separate prediction curves for primiparous
# and multiparous women.
base_primi  <- mf_fit[mf_fit$G1 == "Primip", , drop = FALSE]
base_multi  <- mf_fit[mf_fit$G1 == "Multip", , drop = FALSE]

# 10.4) Build design matrices aligned to fitted fixed-effect columns
# This ensures that predictions align perfectly with the model's fixed effects
build_X <- function(model, newdata){
  mm <- model.matrix(terms(model), newdata)
  Xc <- colnames(getME(model, "X"))
  miss <- setdiff(Xc, colnames(mm))
  if (length(miss))
    mm <- cbind(mm, matrix(0, nrow(mm), length(miss), dimnames = list(NULL, miss)))
  mm[, Xc, drop = FALSE]
}

# Generate a list of design matrices for each ANC value in the grid.
# For each step in the ANC grid we replace the standardized ANC value
# and rebuild the design matrix
make_X_list <- function(base_df, zgrid, model){
  lapply(zgrid, function(z){
    nd <- base_df
    nd$NOC_MX_max <- z
    build_X(model, nd)
  })
}

X_primi <- make_X_list(base_primi, anc_grid_z, logit_mod_gated)
X_multi <- make_X_list(base_multi, anc_grid_z, logit_mod_gated)

# 10.5) Integrate over u ~ N(0, sigma^2) to get population-averaged probabilities
# Compute population-averaged probabilities by integrating over the
# distribution of the random intercept (u ~ Normal(0, sigma²)).
# This produces predictions averaged across possible site-level effects
marginalize_over_u <- function(eta_vec, sd_u, K = 200L){
  if (sd_u <= 0) return(mean(plogis(eta_vec)))
  u <- rnorm(K, 0, sd_u)
  p <- plogis(matrix(eta_vec, nrow = length(eta_vec), ncol = K) +
                rep(u, each = length(eta_vec)))
  mean(p)
}

# Extract the random-intercept standard deviation and fixed-effect coefficients
# from the fitted GLMM
sigma_hat <- as.numeric(attr(VarCorr(logit_mod_gated)[[1]], "stddev"))
beta_hat  <- fixef(logit_mod_gated)

# Generate the predicted probability curve for a given set of design matrices.
# Each ANC value produces a vector of linear predictors that are averaged
# after integrating over the random intercept distribution.
point_curve_from_eta <- function(X_list, beta, sd_u, K = 400L){
  vapply(X_list, function(Xk){
    eta <- as.vector(Xk %*% beta)
    marginalize_over_u(eta, sd_u, K = K)
  }, numeric(1))
}

fig_primi <- point_curve_from_eta(X_primi, beta_hat, sigma_hat)
fig_multi <- point_curve_from_eta(X_multi, beta_hat, sigma_hat)


# 10.6) Uncertainty via MVN draws of beta; attempts profile-based sigma uncertainty (this can be slow)
# Estimate uncertainty by drawing new coefficient vectors from the
# multivariate normal distribution defined by the model estimates
# and their variance-covariance matrix
set.seed(123)
B <- 600L
Vb <- as.matrix(vcov(logit_mod_gated))
beta_sims <- rmvnorm(B, beta_hat, Vb)


# Generate simulated draws of the random-intercept standard deviation.
# If profile confidence intervals are available they are used to approximate
# a log-normal distribution for sigma
get_sigma_draws <- function(mod, B, sigma_hat){
  sigma_sims <- rep(sigma_hat, B)
  ok <- FALSE
  suppressWarnings({
    p <- try(confint(mod, parm = "theta_", method = "profile", oldNames = FALSE), silent = TRUE)
  })
  if (!inherits(p, "try-error")) {
    srow <- grep("^sd_\\(Intercept\\)\\|", rownames(p), value = TRUE)
    if (length(srow) == 1L) {
      lo <- p[srow, 1]; hi <- p[srow, 2]
      if (is.finite(lo) && lo > 0 && is.finite(hi) && hi > 0 && hi > lo) {
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


# For each simulated parameter draw, recompute predicted probabilities
# across the ANC grid. This propagates uncertainty from the model
# coefficients into the predicted curves
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

sim_primi <- mean_probs_from_draws(X_primi, beta_sims, sigma_sims, K_u = 200L)
sim_multi <- mean_probs_from_draws(X_multi, beta_sims, sigma_sims, K_u = 200L)


# Compute 95% percentile confidence intervals across the simulation draws.
ci2 <- function(M)
  rbind(lwr = apply(M, 2, quantile, 0.025),
        upr = apply(M, 2, quantile, 0.975))

pred_fig1b <- rbind(
  data.frame(ANC = anc_grid_raw, Group = "Primiparous",
             Pred = fig_primi, t(ci2(sim_primi))),
  data.frame(ANC = anc_grid_raw, Group = "Multiparous",
             Pred = fig_multi, t(ci2(sim_multi)))
)

# 10.7) Plot Figure 1 and save it
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

fig1 <- plot_curve(pred_fig1b, "", "")
print(fig1)

#ggsave("output/figures/fig1.png", fig1, width = 6, height = 4, dpi = 300)


# 10.8) Brief ANC sanity checks
## does the plot seem reasonable?  
if (isTRUE(attr(sigma_sims, "profile_based"))) {
  message("CIs reflect uncertainty in beta and sigma (sigma via profile-based log-normal draws).")
} else {
  message("CIs reflect uncertainty in beta; sigma held at MLE (no profile CI available).")
}

cat("\nObserved RAW ANC range (model rows):\n"); print(range(anc_raw_used, na.rm = TRUE))
cat("\nCrude Yes/Total at ANC = 1, 4, 8 (overall, model rows):\n")
for (k in c(1,4,8)) {
  ok <- POut$NOC_MX_max[row_idx] == k
  yes <- sum(dat2$home_birthYN[row_idx][ok] == "Yes", na.rm = TRUE)
  tot <- sum(ok, na.rm = TRUE)
  cat(sprintf("ANC=%d: %d/%d (%.1f%%)\n", k, yes, tot, 100*yes/pmax(tot,1)))
}