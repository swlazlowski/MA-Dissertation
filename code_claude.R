###############################################################################
## FULL REPLICATION SCRIPT
## Thesis: "Asymmetric Price Transmission in the UK Petrol Market"
## Original software: Microfit 4.0 | Replication: R
##
## Data: 239 monthly observations, January 1982 – December 2001
## LP = log(net retail price of 4-Star petrol, pence/litre, net of VAT & duty)
## LC = log(Brent crude oil price, USD/barrel, FOB Rotterdam)
## LEX = log(USD/GBP exchange rate)
##
## The CSV has 239 rows (header + 239 data rows → obs 1 = Feb 1982, obs 239 = Dec 2001
## because the thesis says "239 obs from 1982M2 to 2001M12").
## Levels regressions use all 239 rows; ECM/differenced models lose the first obs.
##
## Usage: source("thesis_replication.R")
## Required packages: listed in Section 0. All installed automatically if missing.
################################################################################

# ── 0. PACKAGES ────────────────────────────────────────────────────────────────
pkgs <- c("zoo", "tseries", "lmtest", "sandwich", "car", "strucchange",
          "dynlm", "urca", "forecast", "FinTS", "moments", "MASS")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ── 1. DATA LOAD & TIME SERIES SETUP ──────────────────────────────────────────
# Adjust path if needed
data <- read.csv("data.csv")
# The thesis uses 239 obs from 1982M2 to 2001M12
# (levels: 239 obs; first-differenced: 238 obs from 1982M3)
stopifnot(nrow(data) == 240)

# Create monthly time series (start = February 1982)
LP <- ts(data$LP, start = c(1982, 2), frequency = 12)
LC <- ts(data$LC, start = c(1982, 2), frequency = 12)
LEX <- ts(data$LEX, start = c(1982, 2), frequency = 12)

# Convenience: combined crude-in-GBP = LC + LEX
LCEX <- LC + LEX

# First differences
DLP <- diff(LP)
DLC <- diff(LC)
DLEX <- diff(LEX)
DLCEX <- diff(LCEX)

# Positive / negative first differences (Wolfram / Houck decomposition)
pos <- function(x) pmax(x, 0)
neg <- function(x) pmin(x, 0)

DLP_pos <- pos(DLP); DLP_neg <- neg(DLP)
DLC_pos <- pos(DLC); DLC_neg <- neg(DLC)
DLEX_pos <- pos(DLEX); DLEX_neg <- neg(DLEX)
DLCEX_pos <- pos(DLCEX); DLCEX_neg <- neg(DLCEX)

# Trend variable (1 to 239, aligned to levels)
TREND <- ts(1:240, start = c(1982, 2), frequency = 12)

# Helper: lag a ts by k periods (returns ts of same tsp minus k obs at start)
Lag <- function(x, k = 1) stats::lag(x, -k)

cat("\n====================================================================\n")
cat(" DATA LOADED: ", nrow(data), "monthly obs, Feb 1982 – Dec 2001\n")
cat("====================================================================\n\n")


################################################################################
## TABLE 4 (p.37-38)
## Unit-root tests: DF / ADF for levels and first differences
## Reported in the thesis for LP, LC, LEX and their first differences
################################################################################
cat("── TABLE 4: Unit-root tests ─────────────────────────────────────────\n")

adf_summary <- function(series, name, max_lag = 12){
  # CRDW: run OLS y ~ 1 and compute DW
  m0 <- lm(series ~ 1)
  dw_val <- as.numeric(lmtest::dwtest(m0)$statistic)
  cat(sprintf("\n%s CRDW = %.4f\n", name, dw_val))
  # ADF with lag 0, 6, 12 and IC-chosen lag (BIC)
  for (lag in c(0, 6, 12)) {
    t <- ur.df(series, type = "drift", lags = lag)
    cat(sprintf(" ADF(lag=%2d, drift): tau = %7.4f\n", lag, t@teststat[1]))
    t2 <- ur.df(series, type = "trend", lags = lag)
    cat(sprintf(" ADF(lag=%2d, trend): tau = %7.4f\n", lag, t2@teststat[1]))
  }
  # IC-based selection (BIC up to max_lag)
  bic_vals <- sapply(0:max_lag, function(l) {
    tryCatch(ur.df(series, type = "drift", lags = l)@testreg$coefficients, error = function(e) NULL)
    tryCatch({
      AIC(lm(diff(series) ~ Lag(series, 1) + if (l > 0)
        do.call(cbind, lapply(1:l, function(i) Lag(diff(series), i))) else NULL), k = log(length(series)))
    }, error = function(e) Inf)
  })
  best_lag <- which.min(bic_vals) - 1
  t_best <- ur.df(series, type = "drift", lags = best_lag)
  cat(sprintf(" ADF(IC-chosen lag=%d, drift): tau = %.4f\n", best_lag, t_best@teststat[1]))
}

adf_summary(LP, "LP (log retail price)")
adf_summary(LC, "LC (log crude oil price)")
adf_summary(LEX, "LEX (log exchange rate)")
adf_summary(DLP, "DLP (Δlog retail price)")
adf_summary(DLC, "DLC (Δlog crude oil price)")
adf_summary(DLEX, "DLEX (Δlog exchange rate)")

# Critical values at 5% (MacKinnon 1991): drift ≈ -2.87, trend ≈ -3.43
cat("\nNote: CV at 5%: drift ≈ -2.874, trend ≈ -3.431 (MacKinnon 1991)\n")


################################################################################
## EQUATIONS (16) & (17) (p.39-40)
## Long-run OLS level equations
## (16) LP = β1 + β2*TREND + β3*LC + β4*LEX + ε [with trend]
## (17) LP = β1 + β2*LC + β3*LEX + ε [no trend]
################################################################################
cat("\n\n── EQ (16) & (17): Long-run level equations ──────────────────────────\n")

# --- Equation (16): with trend ---
eq16 <- lm(LP ~ TREND + LC + LEX)
cat("\n=== Equation (16): LP ~ TREND + LC + LEX ===\n")
print(summary(eq16))

# Diagnostic tests (Lagrange Multiplier)
cat(" Serial correlation LM(12):\n")
print(bgtest(eq16, order = 12))
cat(" RESET test:\n")
print(resettest(eq16, power = 2:3, type = "fitted"))
cat(" ARCH(12):\n")
print(FinTS::ArchTest(resid(eq16), lags = 12))
cat(" Heteroscedasticity (BP):\n")
print(bptest(eq16))
cat(" JB normality:\n")
jb16 <- jarque.test(resid(eq16))
print(jb16)
cat(" DW:\n")
print(dwtest(eq16))
cat(" AIC:", AIC(eq16), " SBC (BIC):", BIC(eq16), "\n")

# --- Equation (17): without trend ---
eq17 <- lm(LP ~ LC + LEX)
cat("\n=== Equation (17): LP ~ LC + LEX ===\n")
print(summary(eq17))
print(bgtest(eq17, order = 12))
print(bptest(eq17))
print(jarque.test(resid(eq17)))
cat(" AIC:", AIC(eq17), " BIC:", BIC(eq17), "\n")


################################################################################
## COINTEGRATION TESTS (pp.40-43)
## (a) ADF on residuals from (16) → EG two-step
## (b) CRDW on residuals
## (c) Restricted ECM: t-stat on π (eq 18)
## (d) Unrestricted ECM: Pesaran F-test (eq 19)
## (e) CUSUM and CUSUMSQ stability tests
################################################################################
cat("\n\n── COINTEGRATION TESTS ──────────────────────────────────────────────\n")

resid16 <- ts(resid(eq16), start = c(1982, 2), frequency = 12)

## (a) EG: ADF on residuals (IC-chosen lag = 0 for DF)
cat("\n--- (a) ADF on residuals from eq(16) ---\n")
eg_test <- ur.df(resid16, type = "none", lags = 0)
cat("DF test stat on residuals:", eg_test@teststat, "\n")
cat("Critical values (MacKinnon for 3 regressors+trend):\n")
cat(" 1%: -4.746, 5%: -4.170 (response surfaces, 240 obs)\n")
# AIC/BIC selection on residuals
for (l in 0:3) {
  t <- ur.df(resid16, type = "none", lags = l)
  cat(sprintf(" ADF(lag=%d): tau = %.4f\n", l, t@teststat[1]))
}

## (b) CRDW
cat("\n--- (b) CRDW on residuals ---\n")
dw_resid <- dwtest(lm(resid16 ~ 1))
cat("CRDW:", as.numeric(dw_resid$statistic), " (critical value for 4 regressors ≈ 0.386)\n")

## (c) Restricted ECM eq(18): ΔLP = α + β1*ΔLC + β2*ΔLEX + π*resid(-1) + ε
# Align series (all 238 obs from 1982M3)
resid16_lag <- Lag(resid16, 1)
# Build data frame for ECM (238 obs, 1982M3–2001M12)
ecm_data <- ts.intersect(DLP, DLC, DLEX, resid16_lag)
eq18 <- lm(DLP ~ DLC + DLEX + resid16_lag, data = as.data.frame(ecm_data))
cat("\n=== Equation (18): Restricted ECM ===\n")
print(summary(eq18))
cat(" π (ECT coeff):", coef(eq18)["resid16_lag"], "\n")
cat(" t-stat on ECT:", summary(eq18)$coefficients["resid16_lag", "t value"], "\n")
cat(" CV: 5% = -4.170, 1% = -4.746 (McKinnon 1991)\n")
print(bgtest(eq18, order = 12))
print(bptest(eq18))
print(jarque.test(resid(eq18)))
cat(" AIC:", AIC(eq18), " BIC:", BIC(eq18), "\n")

## (d) Unrestricted ECM eq(19):
## ΔLP = α + β1*ΔLC + β2*ΔLEX + π1*LP(-1) + π2*LC(-1) + π3*LEX(-1) + ε
LP_lag1 <- Lag(LP, 1)
LC_lag1 <- Lag(LC, 1)
LEX_lag1 <- Lag(LEX, 1)
ecm_data2 <- ts.intersect(DLP, DLC, DLEX, LP_lag1, LC_lag1, LEX_lag1)
eq19 <- lm(DLP ~ DLC + DLEX + LP_lag1 + LC_lag1 + LEX_lag1,
           data = as.data.frame(ecm_data2))
cat("\n=== Equation (19): Unrestricted ECM (Pesaran bounds test) ===\n")
print(summary(eq19))
# Pesaran F-test: joint significance of LP(-1), LC(-1), LEX(-1)
pesaran_f <- linearHypothesis(eq19,
                              c("LP_lag1 = 0", "LC_lag1 = 0", "LEX_lag1 = 0"),
                              test = "F")
cat("Pesaran F-test (joint lagged levels):\n")
print(pesaran_f)
cat("CV: 5% upper = 5.972 (k=3, intercept+trend, Pesaran et al. 1996)\n")
cat(" AIC:", AIC(eq19), " BIC:", BIC(eq19), "\n")

## (e) CUSUM and CUSUMSQ
cat("\n--- CUSUM / CUSUMSQ stability tests ---\n")
eq16_struc <- strucchange::efp(LP ~ TREND + LC + LEX, type = "OLS-CUSUM")
eq16_strucSQ <- strucchange::efp(LP ~ TREND + LC + LEX, type = "OLS-MOSUM")
cat("CUSUM: boundary crossings indicate instability\n")
print(sctest(eq16_struc))
# One-step Chow test
ochow <- strucchange::efp(LP ~ TREND + LC + LEX, type = "ME")
cat("One-step Chow (ME) test:\n")
print(sctest(ochow))

# Plots
op <- par(mfrow = c(1, 2))
plot(eq16_struc, main = "CUSUM test (eq 16)")
plot(strucchange::efp(LP ~ TREND + LC + LEX, type = "OLS-CUSUM SQ"),
     main = "CUSUMSQ test (eq 16)")
par(op)


################################################################################
## EQUATION (25) (p.52)
## Long-run equation WITH 5 Pricewatch impulse dummies (Jan–May 1996)
################################################################################
cat("\n\n── EQ (25): Long-run with Pricewatch dummies ────────────────────────\n")

# 5 impulse dummies for Jan–May 1996 (obs 169–173 in the series starting Feb 1982)
n <- length(LP)
dates <- seq(as.Date("1982-02-01"), by = "month", length.out = n)

D1996 <- matrix(0, nrow = n, ncol = 6)
pw_months <- which(format(dates, "%Y-%m") %in%
                     c("1996-01", "1996-02", "1996-03", "1996-04", "1996-05", "1996-06"))
for (i in seq_along(pw_months)) D1996[pw_months[i], i] <- 1
colnames(D1996) <- paste0("D", 1:6)

D1996_ts <- ts(D1996, start = c(1982, 2), frequency = 12)

eq25 <- lm(LP ~ TREND + LC + LEX + D1996_ts)
cat("\n=== Equation (25): LP ~ TREND + LC + LEX +65 dummies ===\n")
print(summary(eq25))
print(bgtest(eq25, order = 12))
print(bptest(eq25))
cat("JB:"); print(jarque.test(resid(eq25)))
cat(" AIC:", AIC(eq25), " BIC:", BIC(eq25), "\n")

# Save residuals for TAR/M-TAR
resid25 <- ts(resid(eq25), start = c(1982, 2), frequency = 12)


################################################################################
## SECTION 3.1 EQUATION (22) (p.46)
## Asymmetric ECM – short-run: Houck/Wolfram approach
## ΔLP = α0 + α1*TREND + β+1,0*ΔLC+ + β1,0*(ΔLC+ + ΔLC(-1)) +
## β+2,0*ΔLEX+ + β2,1*ΔLEX(-1) +
## π1*LP(-1) + π2*(LEX(-1)+LC(-1)) + ε
##
## Thesis reports final parsimonious model after general-to-specific
################################################################################
cat("\n\n── EQ (22): Short-run asymmetric ECM ────────────────────────────────\n")

# Following thesis eq(22) closely:
# ΔLP = 0.23209 + 0.0002869*TREND + 0.16086*ΔLC+_t
# + 0.16357*(ΔLC_t + ΔLC_{t-1})
# + 0.42302*ΔLEX+_{t-1} + 0.34391*ΔLEX_{t-1}
# - 0.15244*LP_{t-1} + 0.06731*(LEX_{t-1}+LC_{t-1})

TREND_d <- Lag(TREND, 0) # aligned to differenced (238 obs)
DLC_pos_L0 <- DLC_pos
DLC_L0 <- DLC
DLC_L1 <- Lag(DLC, 1)
DLEX_pos_L1 <- Lag(DLEX_pos, 1)
DLEX_L1 <- Lag(DLEX, 1)
LP_L1 <- Lag(LP, 1)
LCEX_L1 <- Lag(LCEX, 1) # LC(-1) + LEX(-1)

# Combine contemporaneous and one-lag DLC (thesis constraint: equal coefficients)
DLC_combined <- DLC + Lag(DLC, 1) # ΔLC_t + ΔLC_{t-1}

eq22_data <- ts.intersect(DLP, Lag(TREND, 1), DLC_pos, DLC_combined,
                          Lag(DLEX_pos, 1), Lag(DLEX, 1),
                          Lag(LP, 1), Lag(LCEX, 1))
colnames(eq22_data) <- c("DLP", "TREND", "DLCpos", "DLCcomb",
                         "DLEXpos_L1", "DLEX_L1", "LP_L1", "LCEX_L1")
df22 <- as.data.frame(eq22_data)
eq22 <- lm(DLP ~ TREND + DLCpos + DLCcomb + DLEXpos_L1 + DLEX_L1 + LP_L1 + LCEX_L1,
           data = df22)
cat("\n=== Equation (22): Asymmetric ECM (short-run) ===\n")
print(summary(eq22))
cat(" Serial LM(12):"); print(bgtest(eq22, order = 12))
cat(" RESET:"); print(resettest(eq22, power = 2))
cat(" ARCH(12):"); print(FinTS::ArchTest(resid(eq22), lags = 12))
cat(" BP hetero:"); print(bptest(eq22))
cat(" JB:"); print(jarque.test(resid(eq22)))
cat(" DW:"); print(dwtest(eq22))
cat(" AIC:", AIC(eq22), " BIC:", BIC(eq22), "\n")

# Recursive coefficients on ΔLC+ and ΔLEX+ (Figure 13 & 14)
cat("\nRecursive estimates of asymmetry coefficients (for Figures 13 & 14):\n")
rec <- strucchange::recresid(eq22)
cat(" Recursive residuals computed (", length(rec), "obs)\n")
# Rolling recursive coefficients
rc_full <- data.frame(
  time = time(DLP)[-(1:ncol(model.matrix(eq22)))],
  rec_resid = rec
)


################################################################################
## SECTION 3.2.1 EQUATIONS (24), (26), (28) (pp.50-54)
## ECM with Wolfram segmentation of error-correction term
## Granger & Lee (1989) approach
################################################################################
cat("\n\n── EQ (24), (26), (28): ECM with split ECT ──────────────────────────\n")

## --- Equation (24): simple split ECT, no dummies ---
resid16_lag1 <- Lag(resid16, 1)
ECT_pos_16 <- ts(pmax(lag(resid16, -1), 0), start = c(1982, 3), frequency = 12)
ECT_neg_16 <- ts(pmin(lag(resid16, -1), 0), start = c(1982, 3), frequency = 12)

ecm24_data <- ts.intersect(DLP, DLC, DLEX, ECT_pos_16, ECT_neg_16)
eq24 <- lm(DLP ~ DLC + DLEX + ECT_pos_16 + ECT_neg_16,
           data = as.data.frame(ecm24_data))
cat("\n=== Equation (24): ECM with split ECT (no dummies) ===\n")
print(summary(eq24))
cat(" F-test symmetry (π+ = π-):\n")
print(linearHypothesis(eq24, "ECT_pos_16 = ECT_neg_16"))
cat(" AIC:", AIC(eq24), " BIC:", BIC(eq24), "\n")

## --- Equation (26): split ECT using residuals from eq(25) (with dummies) ---
resid25_lag1_pos <- ts(pmax(lag(resid25, -1), 0), start = c(1982, 3), frequency = 12)
resid25_lag1_neg <- ts(pmin(lag(resid25, -1), 0), start = c(1982, 3), frequency = 12)

ecm26_data <- ts.intersect(DLP, DLC, DLEX, resid25_lag1_pos, resid25_lag1_neg)
eq26 <- lm(DLP ~ DLC + DLEX + resid25_lag1_pos + resid25_lag1_neg,
           data = as.data.frame(ecm26_data))
cat("\n=== Equation (26): ECM with split ECT (dummied residuals) ===\n")
print(summary(eq26))
f26 <- linearHypothesis(eq26, "resid25_lag1_pos = resid25_lag1_neg")
cat(" F-test symmetry (π+ = π-):\n"); print(f26)
# Cook et al. (1998): both must be significant
cat(" t(π+):", summary(eq26)$coef["resid25_lag1_pos", "t value"],
    " t(π-):", summary(eq26)$coef["resid25_lag1_neg", "t value"], "\n")
cat(" AIC:", AIC(eq26), " BIC:", BIC(eq26), "\n")

## --- Equation (28): enriched ECM after general-to-specific ---
# Adding lags of ΔLC, ΔLEX, and split ECT; thesis selects:
# DLC, DLC(-1), DLEX(-1) plus split ECT from resid25
DLC_L1_ts <- Lag(DLC, 1)
DLEX_L1_ts <- Lag(DLEX, 1)
ecm28_data <- ts.intersect(DLP, DLC, DLC_L1_ts, DLEX_L1_ts,
                           resid25_lag1_pos, resid25_lag1_neg)
eq28 <- lm(DLP ~ DLC + DLC_L1_ts + DLEX_L1_ts +
             resid25_lag1_pos + resid25_lag1_neg,
           data = as.data.frame(ecm28_data))
cat("\n=== Equation (28): Enriched ECM (general-to-specific) ===\n")
print(summary(eq28))
f28 <- linearHypothesis(eq28, "resid25_lag1_pos = resid25_lag1_neg")
cat(" F-test symmetry:\n"); print(f28)
cat(" AIC:", AIC(eq28), " BIC:", BIC(eq28), "\n")


################################################################################
## SECTION 3.2.2 TAR AND M-TAR MODELS (pp.57-64)
## Using LCEX = LC + LEX as combined upstream price (in GBP)
##
## Step 1: Long-run equation (29) with LCEX
## Step 2: Eq(32): TAR with threshold = 0
## Step 3: Eq(33): Consistent TAR (Chan's approach, threshold search)
## Step 4: Eq(34): M-TAR with threshold = 0
## Step 5: Eq(35): Consistent M-TAR (threshold search)
################################################################################
cat("\n\n── SECTION 3.2.2: TAR and M-TAR Models ─────────────────────────────\n")

## --- Step 1: Long-run equation (29): LP ~ TREND + LCEX + dummies ---
eq29 <- lm(LP ~ TREND + LCEX + D1996_ts)
cat("\n=== Equation (29): LP ~ TREND + LCEX + dummies ===\n")
print(summary(eq29))
print(bgtest(eq29, order = 12))
print(bptest(eq29))
cat(" JB:"); print(jarque.test(resid(eq29)))
cat(" AIC:", AIC(eq29), " BIC:", BIC(eq29), "\n")

# Test H0: β(LC) = β(LEX) in eq(25) – Wald test for using LCEX
cat("\nWald test: H0: β(LC) = β(LEX) in eq(25) [justifies LCEX]\n")
wtest <- linearHypothesis(eq25, "LC = LEX")
print(wtest)

# Residuals for TAR/M-TAR
resid29 <- ts(resid(eq29), start = c(1982, 2), frequency = 12)

## --- EG test on resid29 ---
cat("\n--- EG cointegration test on resid(29) ---\n")
for (l in 0:2) {
  t <- ur.df(resid29, type = "none", lags = l)
  cat(sprintf(" ADF(lag=%d): tau = %.4f (CV 5%%: -3.823, 1%%: -4.392)\n",
              l, t@teststat[1]))
}

## --- Restricted ECM (30) for LCEX ---
resid29_lag1 <- Lag(resid29, 1)
ecm30_data <- ts.intersect(DLP, DLCEX, resid29_lag1)
eq30 <- lm(DLP ~ DLCEX + resid29_lag1, data = as.data.frame(ecm30_data))
cat("\n=== Equation (30): Restricted ECM (LCEX) ===\n")
print(summary(eq30))
cat(" t on ECT:", summary(eq30)$coef["resid29_lag1","t value"],
    " (CV 5%: -3.820, 1%: -4.392)\n")

## --- Unrestricted ECM (31) ---
LCEX_lag1 <- Lag(LCEX, 1)
LP_lag1_ts <- Lag(LP, 1)
ecm31_data <- ts.intersect(DLP, DLCEX, LP_lag1_ts, LCEX_lag1)
eq31 <- lm(DLP ~ DLCEX + LP_lag1_ts + LCEX_lag1, data = as.data.frame(ecm31_data))
cat("\n=== Equation (31): Unrestricted ECM (LCEX) ===\n")
print(summary(eq31))
f31 <- linearHypothesis(eq31, c("LP_lag1_ts = 0", "LCEX_lag1 = 0"), test = "F")
cat(" Pesaran F (joint lagged levels):\n"); print(f31)
cat(" CV: 5% = 7.423 (k=1, intercept+trend, Pesaran et al. 1996)\n")

################################################################################
## TAR FUNCTIONS
## TAR model: Δu_t = π1 * I(u_{t-1} > τ) * u_{t-1}
## + π2 * (1 - I(u_{t-1} > τ)) * u_{t-1} + ε_t
## M-TAR: threshold on Δu_{t-1} instead of u_{t-1}
################################################################################


################################################################################
##  SECTION 3.2.2  TAR AND M-TAR MODELS  (pp.57-64)
##  Using LCEX = LC + LEX as combined upstream price (in GBP)
##
##  Step 1: Long-run equation (29) with LCEX
##  Step 2: Eq(32): TAR with threshold = 0
##  Step 3: Eq(33): Consistent TAR (Chan's approach, threshold search)
##  Step 4: Eq(34): M-TAR with threshold = 0
##  Step 5: Eq(35): Consistent M-TAR (threshold search)
################################################################################
cat("\n\n── SECTION 3.2.2: TAR and M-TAR Models ─────────────────────────────\n")

## --- Step 1: Long-run equation (29): LP ~ TREND + LCEX + dummies ---
eq29 <- lm(LP ~ TREND + LCEX + D1996_ts)
cat("\n=== Equation (29): LP ~ TREND + LCEX + dummies ===\n")
print(summary(eq29))
print(bgtest(eq29, order = 12))
print(bptest(eq29))
cat("  JB:"); print(jarque.test(resid(eq29)))
cat("  AIC:", AIC(eq29), "  BIC:", BIC(eq29), "\n")

# Test H0: β(LC) = β(LEX) in eq(25) – Wald test for using LCEX
cat("\nWald test: H0: β(LC) = β(LEX) in eq(25) [justifies LCEX]\n")
wtest <- linearHypothesis(eq25, "LC = LEX")
print(wtest)

# Residuals for TAR/M-TAR
resid29 <- ts(resid(eq29), start = c(1982, 2), frequency = 12)

## --- EG test on resid29 ---
cat("\n--- EG cointegration test on resid(29) ---\n")
for (l in 0:2) {
  t <- ur.df(resid29, type = "none", lags = l)
  cat(sprintf("  ADF(lag=%d): tau = %.4f  (CV 5%%: -3.823, 1%%: -4.392)\n",
              l, t@teststat[1]))
}

## --- Restricted ECM (30) for LCEX ---
resid29_lag1 <- Lag(resid29, 1)
ecm30_data <- ts.intersect(DLP, DLCEX, resid29_lag1)
eq30 <- lm(DLP ~ DLCEX + resid29_lag1, data = as.data.frame(ecm30_data))
cat("\n=== Equation (30): Restricted ECM (LCEX) ===\n")
print(summary(eq30))
cat("  t on ECT:", summary(eq30)$coef["resid29_lag1","t value"],
    "  (CV 5%: -3.820, 1%: -4.392)\n")

## --- Unrestricted ECM (31) ---
LCEX_lag1 <- Lag(LCEX, 1)
LP_lag1_ts <- Lag(LP, 1)
ecm31_data <- ts.intersect(DLP, DLCEX, LP_lag1_ts, LCEX_lag1)
eq31 <- lm(DLP ~ DLCEX + LP_lag1_ts + LCEX_lag1, data = as.data.frame(ecm31_data))
cat("\n=== Equation (31): Unrestricted ECM (LCEX) ===\n")
print(summary(eq31))
f31 <- linearHypothesis(eq31, c("LP_lag1_ts = 0", "LCEX_lag1 = 0"), test = "F")
cat("  Pesaran F (joint lagged levels):\n"); print(f31)
cat("  CV: 5% = 7.423 (k=1, intercept+trend, Pesaran et al. 1996)\n")

################################################################################
##  TAR FUNCTIONS
##  TAR model:  Δu_t = π1 * I(u_{t-1} > τ) * u_{t-1}
##                   + π2 * (1 - I(u_{t-1} > τ)) * u_{t-1} + ε_t
##  M-TAR:      threshold on Δu_{t-1} instead of u_{t-1}
################################################################################

run_tar <- function(resids, tau = 0, type = "TAR", label = "") {
  u     <- as.numeric(resids)
  n     <- length(u)
  u_lag <- u[-n]
  du    <- diff(u)
  if (type == "TAR") {
    threshold_var <- u_lag
  } else {  # M-TAR
    threshold_var <- c(NA, diff(u_lag))
    threshold_var <- threshold_var[-1]
    u_lag <- u_lag[-1]
    du    <- du[-1]
  }
  I_ind   <- as.numeric(threshold_var > tau)
  u_pos   <- I_ind * u_lag
  u_neg   <- (1 - I_ind) * u_lag
  fit     <- lm(du ~ 0 + u_pos + u_neg)
  pi1     <- coef(fit)["u_pos"]
  pi2     <- coef(fit)["u_neg"]
  tstat   <- summary(fit)$coef[, "t value"]
  t_max   <- max(tstat)         # t-max statistic
  fstat_phi <- summary(fit)$fstatistic[1]  # ϕ: joint H0: π1=π2=0
  f_sym   <- linearHypothesis(fit, "u_pos = u_neg")$F[2]  # H0: π1=π2
  rss     <- sum(resid(fit)^2)
  cat(sprintf("\n=== %s model (τ = %.5f) ===\n", label, tau))
  print(summary(fit))
  cat(sprintf("  π1 = %.5f  π2 = %.5f\n", pi1, pi2))
  cat(sprintf("  t-max = %.4f  (CV 5%%: -2.12, 1%%: -2.53 for TAR; -1.99, -2.45 for M-TAR)\n",
              t_max))
  cat(sprintf("  ϕ (joint π1=π2=0) = %.4f  (CV 5%%: 5.87, 1%%: 8.04)\n", fstat_phi))
  cat(sprintf("  F(π1=π2) = %.4f  p-value = %.4f\n",
              f_sym, linearHypothesis(fit, "u_pos = u_neg")$Pr[2]))
  cat(sprintf("  RSS = %.5f\n", rss))
  invisible(list(fit = fit, pi1 = pi1, pi2 = pi2, t_max = t_max,
                 phi = fstat_phi, f_sym = f_sym, rss = rss, tau = tau))
}

## --- Chan's threshold search (consistent TAR / M-TAR) ---
chan_search <- function(resids, type = "TAR", trim = 0.15) {
  u   <- as.numeric(resids)
  n   <- length(u)
  u_lag <- u[-n]
  du  <- diff(u)
  if (type == "MTAR") {
    th_var <- c(NA, diff(u_lag)); th_var <- th_var[-1]
    u_lag  <- u_lag[-1]; du <- du[-1]
  } else {
    th_var <- u_lag
  }
  # Sort and trim
  sorted <- sort(th_var)
  lo <- floor(trim * length(sorted)) + 1
  hi <- floor((1 - trim) * length(sorted))
  candidates <- sorted[lo:hi]
  rss_vec <- sapply(candidates, function(tau) {
    I_ind <- as.numeric(th_var > tau)
    u_pos <- I_ind * u_lag; u_neg <- (1 - I_ind) * u_lag
    sum(resid(lm(du ~ 0 + u_pos + u_neg))^2)
  })
  best_tau <- candidates[which.min(rss_vec)]
  best_tau
}

## --- Equation (32): TAR, τ = 0 ---
cat("\n\n── EQ (32): TAR model, threshold = 0 ───────────────────────────────\n")
res32 <- run_tar(resid29, tau = 0, type = "TAR", label = "Eq(32) TAR")

## --- Equation (33): Consistent TAR (Chan's τ) ---
cat("\n── EQ (33): Consistent TAR (Chan's method) ──────────────────────────\n")
tau33 <- chan_search(resid29, type = "TAR")
cat(sprintf("Chan's optimal threshold τ = %.5f\n", tau33))
res33 <- run_tar(resid29, tau = tau33, type = "TAR", label = "Eq(33) Consistent TAR")
cat("Thesis reports τ = -0.05969, RSS = 0.83589 (ordinary TAR RSS = 0.85388)\n")

## --- Equation (34): M-TAR, τ = 0 ---
cat("\n── EQ (34): M-TAR model, threshold = 0 ─────────────────────────────\n")
res34 <- run_tar(resid29, tau = 0, type = "MTAR", label = "Eq(34) M-TAR")

## --- Equation (35): Consistent M-TAR ---
cat("\n── EQ (35): Consistent M-TAR (Chan's method) ───────────────────────\n")
tau35 <- chan_search(resid29, type = "MTAR")
cat(sprintf("Chan's optimal threshold (M-TAR) τ = %.5f\n", tau35))
res35 <- run_tar(resid29, tau = tau35, type = "MTAR", label = "Eq(35) Consistent M-TAR")
cat("Thesis reports τ = -0.0034, RSS = 0.82965\n")


################################################################################
##  SECTION 4  EXPLANATIONS FOR ASYMMETRY  (pp.68-79)
##  Table 5: Asymmetry × volatility of petrol prices
##  Tables 7 & 8: Asymmetry × change in cost share
##  Tables 9 & 10: Asymmetry × volatility of input prices (Peltzman test)
##  Tables 11 & 12: Asymmetry × drift of input prices (menu costs)
################################################################################
cat("\n\n── SECTION 4: Explanations for asymmetry ───────────────────────────\n")

# Build a data frame of aligned levels (all 239 obs)
p_lev   <- as.numeric(LP)    # log retail
c_lev   <- as.numeric(LC)    # log crude USD
ex_lev  <- as.numeric(LEX)   # log exchange rate
cex_lev <- as.numeric(LCEX)  # log crude GBP

# Price ratio: upstream cost share = exp(LC + LEX) / exp(LP) = exp(LCEX - LP)
cost_share <- exp(cex_lev - p_lev)
d_cost_share <- c(NA, diff(cost_share))

# Monthly changes in petrol price
dlp_vec  <- c(NA, diff(p_lev))
dlc_vec  <- c(NA, diff(c_lev))
dlex_vec <- c(NA, diff(ex_lev))
dlcex_vec <- c(NA, diff(cex_lev))

dlc_pos_v  <- pmax(dlc_vec,  0)
dlex_pos_v <- pmax(dlex_vec, 0)

# Align everything as a data frame (238 rows: obs 2..239)
align_df <- function() {
  n239 <- 239
  data.frame(
    dlp    = dlp_vec[2:n239],
    dlc    = dlc_vec[2:n239],
    dlex   = dlex_vec[2:n239],
    dlcex  = dlcex_vec[2:n239],
    dlc_pos  = dlc_pos_v[2:n239],
    dlex_pos = dlex_pos_v[2:n239],
    lp_lag  = p_lev[1:(n239-1)],
    lcex_lag = cex_lev[1:(n239-1)],
    trend   = 1:(n239-1),
    cost_share_chg = d_cost_share[2:n239]
  )
}
df_main <- align_df()

## ── Table 5: Asymmetry & volatility of petrol prices ──────────────────────
cat("\n--- TABLE 5: Asymmetry and volatility of petrol prices ---\n")
sd_dlp <- sd(df_main$dlp, na.rm = TRUE)
thresholds_pct <- c(1.00, 0.50, 0.33, 0.25)

for (pct in thresholds_pct) {
  theta <- pct * sd_dlp
  D_vol <- as.numeric(abs(df_main$dlp) < theta)   # 1 = low volatility period
  # Interacted model: D*ΔLC+ and ΔLC+ (separate effects)
  fit_c <- lm(dlp ~ trend + I(D_vol * dlc_pos) + dlc_pos +
                I(D_vol * dlex_pos) + dlex_pos +
                lp_lag + lcex_lag,
              data = df_main)
  cat(sprintf("\n  Threshold = %.0f%% sd (%.4f):\n", pct*100, theta))
  cat(sprintf("    D*Δ+lnc = %.5f  (t=%.3f)\n",
              coef(fit_c)["I(D_vol * dlc_pos)"],
              summary(fit_c)$coef["I(D_vol * dlc_pos)","t value"]))
  cat(sprintf("    Δ+lnc   = %.5f  (t=%.3f)\n",
              coef(fit_c)["dlc_pos"],
              summary(fit_c)$coef["dlc_pos","t value"]))
  cat(sprintf("    D*Δ+lex = %.5f  (t=%.3f)\n",
              coef(fit_c)["I(D_vol * dlex_pos)"],
              summary(fit_c)$coef["I(D_vol * dlex_pos)","t value"]))
  cat(sprintf("    Δ+lex   = %.5f  (t=%.3f)\n",
              coef(fit_c)["dlex_pos"],
              summary(fit_c)$coef["dlex_pos","t value"]))
  # Wald: H0: sum = 0 (symmetric when volatile)
  w <- tryCatch(linearHypothesis(fit_c,
                                 "I(D_vol * dlc_pos) + dlc_pos = 0"), error = function(e) NULL)
  if (!is.null(w)) cat(sprintf("    Wald(Σ=0): F=%.4f  p=%.4f\n",
                               w$F[2], w$Pr[2]))
}

## ── Tables 7 & 8: Asymmetry × change in cost share ────────────────────────
cat("\n--- TABLES 7 & 8: Asymmetry and change in cost share ---\n")
# k-period moving averages of cost share change
for (k in c(1, 2, 3, 4, 6, 9, 12)) {
  ma_cs <- filter(df_main$cost_share_chg, rep(1/k, k), sides = 1)
  D_cs  <- as.numeric(!is.na(ma_cs) & ma_cs > 0)   # 1 = rising cost share
  fit_cs <- lm(dlp ~ trend + I(D_cs * dlc_pos) + dlc_pos +
                 I(D_cs * dlex_pos) + dlex_pos +
                 lp_lag + lcex_lag, data = df_main)
  cat(sprintf("\n  MA order k=%2d:\n", k))
  cat(sprintf("    D*Δ+lnc = %.5f  Δ+lnc = %.5f\n",
              coef(fit_cs)["I(D_cs * dlc_pos)"],
              coef(fit_cs)["dlc_pos"]))
  # Wald H0: coefs are equal in magnitude, opposite signs (symmetric)
  w2 <- tryCatch(linearHypothesis(fit_cs,
                                  "I(D_cs * dlc_pos) + dlc_pos = 0"), error = function(e) NULL)
  if (!is.null(w2)) cat(sprintf("    Wald(Σ=0): F=%.4f  p=%.4f\n",
                                w2$F[2], w2$Pr[2]))
}

## ── Tables 9 & 10: Asymmetry × input price volatility (Peltzman) ──────────
cat("\n--- TABLES 9 & 10: Asymmetry and input price volatility ---\n")
for (var_name in c("dlc", "dlcex")) {
  cat(sprintf("\n  Input variable: %s\n", var_name))
  inp <- df_main[[var_name]]
  for (k in c(4, 6, 12, 24)) {
    ma_inp <- filter(inp, rep(1/k, k), sides = 1)
    for (xi_pct in c(1.00, 0.50, 0.33, 0.25)) {
      dev   <- inp - ma_inp           # deviation from MA
      xi    <- xi_pct * sd(dev, na.rm = TRUE)
      D_vol <- as.numeric(!is.na(dev) & abs(dev) < xi)  # 1 = low vol
      fit_v <- lm(dlp ~ trend + I(D_vol * dlc_pos) + dlc_pos +
                    I(D_vol * dlex_pos) + dlex_pos +
                    lp_lag + lcex_lag, data = df_main)
      c1 <- coef(fit_v)["I(D_vol * dlc_pos)"]
      c2 <- coef(fit_v)["dlc_pos"]
      cat(sprintf("    k=%2d, xi=%.0f%%sd: D*Δ+=%7.5f  Δ+=%7.5f\n",
                  k, xi_pct*100, c1, c2))
    }
  }
}

## ── Tables 11 & 12: Asymmetry × drift of input prices (menu costs) ─────────
cat("\n--- TABLES 11 & 12: Asymmetry and drift of input prices ---\n")
for (var_name in c("dlc", "dlcex")) {
  cat(sprintf("\n  Drift variable: %s\n", var_name))
  inp <- df_main[[var_name]]
  for (k in c(4, 6, 12, 18)) {
    ma_drift <- filter(inp, rep(1/k, k), sides = 1)
    D_drift  <- as.numeric(!is.na(ma_drift) & ma_drift > 0)  # 1 = upward drift
    fit_d <- lm(dlp ~ trend + I(D_drift * dlc_pos) + dlc_pos +
                  I(D_drift * dlex_pos) + dlex_pos +
                  lp_lag + lcex_lag, data = df_main)
    c1 <- coef(fit_d)["I(D_drift * dlc_pos)"]
    c2 <- coef(fit_d)["dlc_pos"]
    t1 <- summary(fit_d)$coef["I(D_drift * dlc_pos)", "t value"]
    t2 <- summary(fit_d)$coef["dlc_pos", "t value"]
    cat(sprintf("    k=%2d: D*Δ+=%7.5f (t=%.3f)  Δ+=%7.5f (t=%.3f)\n",
                k, c1, t1, c2, t2))
  }
}


################################################################################
##  REPLICATION SUMMARY TABLE
##  Collects key numbers for comparison against thesis & Microfit output
################################################################################
cat("\n\n")
cat("====================================================================\n")
cat(" REPLICATION SUMMARY – compare against thesis / Microfit printouts\n")
cat("====================================================================\n\n")

cat("--- EQ(16) long-run [with trend] ---\n")
cf16 <- coef(eq16)
cat(sprintf("  Intercept=%.4f  Trend=%.6f  LC=%.5f  LEX=%.5f\n",
            cf16[1], cf16[2], cf16[3], cf16[4]))
cat(sprintf("  R2=%.3f  DW=%.2f  AIC=%.2f  BIC=%.2f\n",
            summary(eq16)$r.squared, as.numeric(dwtest(eq16)$statistic),
            AIC(eq16), BIC(eq16)))
cat("  Thesis reports: α=1.3894  trend=0.003729  LC=0.49011  LEX=0.37999\n\n")

cat("--- EQ(17) long-run [no trend] ---\n")
cf17 <- coef(eq17)
cat(sprintf("  Intercept=%.4f  LC=%.5f  LEX=%.5f\n", cf17[1], cf17[2], cf17[3]))
cat("  Thesis reports: α=1.5126  LC=0.46407  LEX=0.37947\n\n")

cat("--- EQ(25) long-run [with dummies] ---\n")
cf25 <- coef(eq25)
cat(sprintf("  Intercept=%.4f  Trend=%.7f  LC=%.5f  LEX=%.5f\n",
            cf25[1], cf25[2], cf25[3], cf25[4]))
cat("  Thesis reports: α=1.4008  trend=0.0004374  LC=0.48936  LEX=0.40411\n\n")

cat("--- EQ(29) long-run [LCEX, with dummies] ---\n")
cf29 <- coef(eq29)
cat(sprintf("  Intercept=%.4f  Trend=%.6f  LCEX=%.5f\n",
            cf29[1], cf29[2], cf29[3]))
cat("  Thesis reports: α=1.4864  trend=0.003600  LCEX=0.47251\n\n")

cat("--- EQ(32) TAR τ=0 ---\n")
cat(sprintf("  π1(pos)=%.5f  π2(neg)=%.5f\n", res32$pi1, res32$pi2))
cat(sprintf("  t-max=%.4f  ϕ=%.4f  F(sym)=%.4f  RSS=%.5f\n",
            res32$t_max, res32$phi, res32$f_sym, res32$rss))
cat("  Thesis reports: π1=-0.21648  π2=-0.32466  t-max=-3.4599  F(sym)=1.4976\n\n")

cat("--- EQ(33) Consistent TAR ---\n")
cat(sprintf("  τ=%.5f  π1=%.5f  π2=%.5f\n", res33$tau, res33$pi1, res33$pi2))
cat(sprintf("  F(sym)=%.4f  RSS=%.5f\n", res33$f_sym, res33$rss))
cat("  Thesis reports: τ=-0.05969  π1=-0.17884  π2=-0.38302  F(sym)=5.366\n\n")

cat("--- EQ(34) M-TAR τ=0 ---\n")
cat(sprintf("  π1=%.5f  π2=%.5f  F(sym)=%.4f\n", res34$pi1, res34$pi2, res34$f_sym))
cat("  Thesis reports: π1=-0.15352  π2=-0.38019  F(sym)=6.68\n\n")

cat("--- EQ(35) Consistent M-TAR ---\n")
cat(sprintf("  τ=%.5f  π1=%.5f  π2=%.5f  F(sym)=%.4f\n",
            res35$tau, res35$pi1, res35$pi2, res35$f_sym))
cat("  Thesis reports: τ=-0.0034  π1=-0.15026  π2=-0.38923  Wald=7.4643\n\n")

cat("====================================================================\n")
cat(" Notes for grading:\n")
cat("  1. Minor numerical differences (±0.001) are normal due to\n")
cat("     Microfit vs R floating-point implementation.\n")
cat("  2. Large discrepancies signal errors in the student's reported output.\n")
cat("  3. The Microfit .doc file ends with:\n")
cat("     'ok  zgadza sie z wydrukami i z estymacja' [Polish: 'ok agrees\n")
cat("     with printouts and estimation'] — confirming own verification.\n")
cat("====================================================================\n")

