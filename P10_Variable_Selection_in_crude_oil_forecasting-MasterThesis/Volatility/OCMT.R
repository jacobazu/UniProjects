

library(dplyr)
library(tidyr)
library(ggplot2)
library(forecast)   
library(rugarch)

# ----------------------------------------------------------------------
# 0. User settings
# ----------------------------------------------------------------------
DATA_FILE <- "VolatilityData.csv"
DATE_VAR  <- "date"
TARGET_VAR <- "LV"

WINDOW_SIZE <- 60
pval_grid   <- c(0.01, 0.05, 0.10)
delta1_grid <- 1
delta2_grid <- 2
CV_SCALE <- 1

TRAIN_END <- as.Date("2010-12-01")
VAL_END   <- as.Date("2013-12-01")
TEST_END  <- as.Date("2024-02-01")

# ----------------------------------------------------------------------
# 1. Helper functions 
# ----------------------------------------------------------------------
cp_threshold <- function(pval, delta, n, cscale = 1) {
  if (n <= 0) stop("n must be positive.")
  qnorm(1 - pval / (2 * cscale * (n ^ delta)))
}

safe_lm_tstat <- function(y, x, Xsel = NULL, Xcond = NULL) {
  y <- as.numeric(y)
  x <- as.numeric(x)
  n <- length(y)
  if (length(x) != n) stop("x and y must have same length.")
  
  if (is.null(Xsel)) Xsel <- matrix(nrow = n, ncol = 0)
  else Xsel <- as.matrix(Xsel)
  if (is.null(Xcond)) Xcond <- matrix(nrow = n, ncol = 0)
  else Xcond <- as.matrix(Xcond)
  
  df <- data.frame(y = y, x = x)
  if (ncol(Xsel) > 0) {
    ctrl1 <- as.data.frame(Xsel)
    names(ctrl1) <- paste0("s", seq_len(ncol(ctrl1)))
    df <- cbind(df, ctrl1)
  }
  if (ncol(Xcond) > 0) {
    ctrl2 <- as.data.frame(Xcond)
    names(ctrl2) <- paste0("z", seq_len(ncol(ctrl2)))
    df <- cbind(df, ctrl2)
  }
  
  fit <- tryCatch(lm(y ~ ., data = df), error = function(e) NULL)
  if (is.null(fit)) return(-Inf)
  sm <- summary(fit)
  if (!("x" %in% rownames(sm$coefficients))) return(-Inf)
  abs(sm$coefficients["x", "t value"])
}

ocmt_select <- function(y, X, Xcond = NULL, pval = 0.05,
                        delta1 = 1, delta2 = 2, cscale = 1,
                        max_stages = 50) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  if (length(y) != n) stop("length(y) must equal nrow(X).")
  if (is.null(Xcond)) Xcond <- matrix(nrow = n, ncol = 0)
  else Xcond <- as.matrix(Xcond)
  
  thr1 <- cp_threshold(pval, delta1, p, cscale)
  thr2 <- cp_threshold(pval, delta2, p, cscale)
  
  selected <- rep(FALSE, p)
  stage <- 1L
  repeat {
    active <- which(!selected)
    if (length(active) == 0) break
    thr <- if (stage == 1L) thr1 else thr2
    tstats <- rep(-Inf, p)
    for (j in active) {
      Xsel <- if (any(selected)) X[, selected, drop = FALSE] else NULL
      tstats[j] <- safe_lm_tstat(y, X[, j], Xsel, Xcond)
    }
    sel_stage <- active[is.finite(tstats[active]) & (tstats[active] > thr)]
    if (length(sel_stage) == 0) break
    selected[sel_stage] <- TRUE
    stage <- stage + 1L
    if (stage > max_stages) break
  }
  list(selected = selected, selected_names = colnames(X)[selected])
}

forecast_ocmt <- function(y_train, X_train, X_test, selected, Xcond_train = NULL, Xcond_test = NULL) {
  if (sum(selected) == 0) return(mean(y_train))
  Xsel_train <- X_train[, selected, drop = FALSE]
  Xsel_test  <- X_test[, selected, drop = FALSE]
  df_train <- data.frame(y = y_train, Xsel_train)
  fit <- tryCatch(lm(y ~ ., data = df_train), error = function(e) NULL)
  if (is.null(fit)) return(mean(y_train))
  newdata <- as.data.frame(as.list(Xsel_test[1, , drop = TRUE]))
  names(newdata) <- colnames(Xsel_train)
  pred <- suppressWarnings(predict(fit, newdata = newdata))[1]
  if (!is.finite(pred)) mean(y_train) else pred
}

dm_test <- function(y, f_model, f_bench) {
  e_model <- y - f_model
  e_bench <- y - f_bench
  d <- e_bench^2 - e_model^2
  n <- length(d)
  stat <- mean(d) / (sd(d) / sqrt(n))
  pval <- 1 - pnorm(stat)
  list(stat = stat, p_value = pval)
}

cw_test <- function(y, f_model, f_bench) {
  e_b <- y - f_bench
  e_m <- y - f_model
  d_cw <- e_b^2 - (e_m^2 - (f_bench - f_model)^2)
  n <- length(d_cw)
  stat <- mean(d_cw) / (sd(d_cw) / sqrt(n))
  pval <- 1 - pnorm(stat)
  list(stat = stat, p_value = pval)
}

# OCMT evaluation 
run_rolling_ocmt <- function(y_all, X_all, dates,
                             start_idx, end_idx,
                             pval, delta1, delta2,
                             window_size = 60,          # rolling window for selection
                             Xcond_all = NULL) {
  n_eval <- end_idx - start_idx + 1
  p <- ncol(X_all)
  f_model <- numeric(n_eval)
  f_bench <- numeric(n_eval)
  y_oos   <- numeric(n_eval)
  sel_mat <- matrix(0L, nrow = n_eval, ncol = p)
  colnames(sel_mat) <- colnames(X_all)
  
  for (ii in seq_len(n_eval)) {
    t <- start_idx + ii - 1
    

    sel_start <- 1
    sel_end <- t - 1
    sel_idx <- sel_start:sel_end
    y_sel <- y_all[sel_idx]
    X_sel <- X_all[sel_idx, , drop = FALSE]
    

    bench_start <- 1
    bench_end <- t - 1
    bench_idx <- bench_start:bench_end
    y_bench <- y_all[bench_idx]
    
    Xcond_sel <- if (!is.null(Xcond_all)) Xcond_all[sel_idx, , drop = FALSE] else NULL
    Xcond_bench <- NULL  # not used for benchmark
    

    sel <- ocmt_select(y_sel, X_sel, Xcond_sel, pval, delta1, delta2, CV_SCALE)
    
    # Forecast using selected predictors 
    f1 <- forecast_ocmt(y_sel, X_sel, X_all[t, , drop = FALSE], sel$selected, Xcond_sel, NULL)
    f0 <- mean(y_bench, na.rm = TRUE)  
    
    y_oos[ii] <- y_all[t]
    f_bench[ii] <- f0
    f_model[ii] <- f1
    sel_mat[ii, ] <- as.integer(sel$selected)
  }
  
  list(y = y_oos, f_model = f_model, f_bench = f_bench,
       selection_matrix = sel_mat)
}

# ----------------------------------------------------------------------
# 2. Load and prepare data
# ----------------------------------------------------------------------
data <- read.csv(DATA_FILE)
data[[DATE_VAR]] <- as.Date(data[[DATE_VAR]])
predictors <- setdiff(names(data), c(DATE_VAR, TARGET_VAR))
keep_cols <- c(DATE_VAR, TARGET_VAR, predictors)
df <- data %>%
  select(all_of(keep_cols)) %>%
  filter(complete.cases(.)) %>%
  arrange(.data[[DATE_VAR]])

dates <- df[[DATE_VAR]]
y_all <- df[[TARGET_VAR]]
X_all <- as.matrix(df[, predictors, drop = FALSE])
colnames(X_all) <- predictors

train_end_idx <- max(which(dates <= TRAIN_END))
val_end_idx   <- max(which(dates <= VAL_END))
test_end_idx  <- max(which(dates <= TEST_END))
if (!is.finite(train_end_idx) || !is.finite(val_end_idx) || !is.finite(test_end_idx))
  stop("Date ranges not found in data.")

cat("Training   :", as.character(min(dates[1:train_end_idx])), "to", as.character(max(dates[1:train_end_idx])), "\n")
cat("Validation :", as.character(dates[train_end_idx+1]), "to", as.character(dates[val_end_idx]), "\n")
cat("Test       :", as.character(dates[val_end_idx+1]), "to", as.character(dates[test_end_idx]), "\n")

# ----------------------------------------------------------------------
# 3. Grid search on validation period
# ----------------------------------------------------------------------
grid <- expand.grid(pval = pval_grid, delta1 = delta1_grid, delta2 = delta2_grid,
                    stringsAsFactors = FALSE)
grid <- grid[grid$delta2 > grid$delta1, ]
grid$rel_rmsfe <- NA_real_

for (g in seq_len(nrow(grid))) {
  cat("Evaluating (pval, d1, d2) =", grid$pval[g], ",", grid$delta1[g], ",", grid$delta2[g], "... ")
  res_val <- run_rolling_ocmt(y_all, X_all, dates,
                              start_idx = train_end_idx + 1,
                              end_idx = val_end_idx,
                              pval = grid$pval[g],
                              delta1 = grid$delta1[g],
                              delta2 = grid$delta2[g],
                              window_size = WINDOW_SIZE)
  rmsfe_model <- sqrt(mean((res_val$y - res_val$f_model)^2))
  rmsfe_bench <- sqrt(mean((res_val$y - res_val$f_bench)^2))
  grid$rel_rmsfe[g] <- rmsfe_model / rmsfe_bench
  cat("rel RMSFE =", round(grid$rel_rmsfe[g], 4), "\n")
}

best <- grid[which.min(grid$rel_rmsfe), ]
cat("\nBest hyperparameters:\n")
print(best)

# ----------------------------------------------------------------------
# 4. Final evaluation on test period
# ----------------------------------------------------------------------
final <- run_rolling_ocmt(y_all, X_all, dates,
                          start_idx = val_end_idx + 1,
                          end_idx = test_end_idx,
                          pval = best$pval,
                          delta1 = best$delta1,
                          delta2 = best$delta2,
                          window_size = WINDOW_SIZE)

actual <- final$y
pred_ocmt <- final$f_model
pred_bench <- final$f_bench

# Point forecast metrics
mspe_ocmt <- mean((actual - pred_ocmt)^2)
mspe_bench <- mean((actual - pred_bench)^2)
oos_r2 <- 1 - mspe_ocmt / mspe_bench
rmsfe_ocmt <- sqrt(mspe_ocmt)
rmsfe_bench <- sqrt(mspe_bench)
rel_rmsfe <- rmsfe_ocmt / rmsfe_bench

# Directional metrics on changes
actual_diff <- diff(actual)
pred_diff <- diff(pred_ocmt)
success <- mean(sign(actual_diff) == sign(pred_diff))
pt <- DACTest(pred_diff, actual_diff, test = "PT")

dm <- dm_test(actual, pred_ocmt, pred_bench)
cw <- cw_test(actual, pred_ocmt, pred_bench)

cat("\n========== FINAL TEST RESULTS ==========\n")
cat("RMSFE (OCMT):", rmsfe_ocmt, "\n")
cat("RMSFE (Benchmark):", rmsfe_bench, "\n")
cat("Relative RMSFE:", rel_rmsfe, "\n")
cat("OOS R²:", oos_r2, "\n")
cat("Success ratio (on changes):", success, "\n")
cat("Diebold-Mariano test: stat =", dm$stat, " p-value =", dm$p_value, "\n")
cat("Pesaran-Timmermann test (on changes):\n")
print(pt)
cat("Clark-West test: stat =", cw$stat, " p-value =", cw$p_value, "\n")

# ----------------------------------------------------------------------
# 5. CSPE difference plot (OCMT vs Historical mean)
# ----------------------------------------------------------------------
cspe_diff <- cumsum((actual - pred_bench)^2 - (actual - pred_ocmt)^2)
cspe_df <- data.frame(Date = dates[(val_end_idx + 1):test_end_idx], CSPE_Diff = cspe_diff)
p_cspe <- ggplot(cspe_df, aes(x = Date, y = CSPE_Diff)) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  theme_classic() +
  labs(title = paste("Cumulative squared prediction error difference (OCMT vs Historical mean)\n",
                     "pval =", best$pval, ", δ1 =", best$delta1, ", δ2 =", best$delta2),
       x = "Time", y = "CSPE (Benchmark - OCMT)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cspe)

# ----------------------------------------------------------------------
# 6. Model size over time
# ----------------------------------------------------------------------
model_size <- rowSums(final$selection_matrix)
test_dates <- dates[(val_end_idx + 1):test_end_idx]
model_size_df <- data.frame(Date = test_dates, ModelSize = model_size)
p_size <- ggplot(model_size_df, aes(x = Date, y = ModelSize)) +
  geom_line(color = "darkred", linewidth = 0.8) +
  geom_point(size = 0.8, color = "darkred") +
  theme_classic() +
  labs(title = "OCMT model size over time",
       x = "Time", y = "Number of selected predictors") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_size)
cat("\nModel size statistics:\n")
cat("  Mean :", round(mean(model_size), 2), "\n")
cat("  Median:", median(model_size), "\n")
cat("  Min  :", min(model_size), "\n")
cat("  Max  :", max(model_size), "\n")

# ----------------------------------------------------------------------
# 7. Variable selection frequencies
# ----------------------------------------------------------------------
freq <- colMeans(final$selection_matrix)
freq_pct <- sort(100 * freq, decreasing = TRUE)
cat("\n--- Variable selection frequencies (%) ---\n")
print(round(freq_pct, 1))

ordered_vars <- names(freq_pct)
if (length(ordered_vars) > 0) {
  sel_plot <- final$selection_matrix[, ordered_vars, drop = FALSE]
  heat_df <- as.data.frame(sel_plot)
  heat_df$time <- test_dates
  long_df <- heat_df %>%
    pivot_longer(cols = all_of(ordered_vars), names_to = "variable", values_to = "selected")
  p_heat <- ggplot(long_df, aes(x = time, y = variable, fill = factor(selected))) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "white", "1" = "blue"),
                      labels = c("Not selected", "Selected"), name = "") +
    labs(title = "OCMT Selection Heatmap", x = "Time", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 7), panel.grid = element_blank())
  print(p_heat)
}

# ----------------------------------------------------------------------
# VaR savings using Filtered Historical Simulation (FHS)
# ----------------------------------------------------------------------
# Load futures returns from Zhang dataset
zhang_data <- read.csv("ZhangData32_3.0.csv")
zhang_data$date <- as.Date(zhang_data$date)
zhang_data <- zhang_data %>%
  arrange(date) %>%
  mutate(futures_ret = (Price_Futures - lag(Price_Futures)) / lag(Price_Futures)) %>%
  select(date, futures_ret)

# Align forecasts with test dates
var_data <- data.frame(date = test_dates,
                       LV_forecast = pred_ocmt,
                       LV_bench = pred_bench) %>%
  left_join(zhang_data, by = "date") %>%
  filter(complete.cases(.))

# Convert LV forecasts to monthly standard deviation
var_data <- var_data %>%
  mutate(
    sigma_model = exp(LV_forecast / 2),
    sigma_bench = exp(LV_bench / 2)
  )

# FHS: rolling empirical quantile of standardized returns
window_size <- 60   # months
empirical_quantile <- rep(NA_real_, nrow(var_data))

for (i in 1:nrow(var_data)) {
  past_returns <- var_data$futures_ret[1:(i-1)]
  if (length(past_returns) >= window_size) {
    recent <- tail(past_returns, window_size)
    std_recent <- (recent - mean(recent, na.rm = TRUE)) / sd(recent, na.rm = TRUE)
    empirical_quantile[i] <- quantile(std_recent, probs = 0.05, na.rm = TRUE)
  } else {
    empirical_quantile[i] <- qnorm(0.05)   # fallback
  }
}

# VaR = - sigma * empirical_quantile
var_data <- var_data %>%
  mutate(
    VaR_model = - sigma_model * empirical_quantile,
    VaR_bench = - sigma_bench * empirical_quantile,
    VaR_saving = VaR_bench - VaR_model
  )

cat("\n--- VaR savings (95% confidence, FHS) ---\n")
cat("Average VaR (benchmark):", round(mean(var_data$VaR_bench, na.rm = TRUE), 6), "\n")
cat("Average VaR (model):    ", round(mean(var_data$VaR_model, na.rm = TRUE), 6), "\n")
cat("Average VaR saving:     ", round(mean(var_data$VaR_saving, na.rm = TRUE), 6), "\n")
pct_saving <- 100 * mean(var_data$VaR_saving, na.rm = TRUE) / mean(var_data$VaR_bench, na.rm = TRUE)
cat("Average VaR saving (%): ", round(pct_saving, 2), "%\n")

# Violation rates
var_data <- var_data %>%
  mutate(
    exceed_model = ifelse(futures_ret < -VaR_model, 1, 0),
    exceed_bench = ifelse(futures_ret < -VaR_bench, 1, 0)
  )
violation_rate_model <- mean(var_data$exceed_model, na.rm = TRUE)
violation_rate_bench <- mean(var_data$exceed_bench, na.rm = TRUE)
cat("\nVaR violation rate (model):   ", round(violation_rate_model, 4), "\n")
cat("VaR violation rate (benchmark):", round(violation_rate_bench, 4), "\n")
cat("Nominal level (alpha):         0.05\n")

# Plot VaR savings
p_var <- ggplot(var_data, aes(x = date, y = VaR_saving)) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  theme_classic() +
  labs(title = "VaR savings (OCMT vs Historical mean) - Filtered Historical Simulation",
       x = "Time", y = "VaR reduction") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_var)

# ----------------------------------------------------------------------
# Straddle strategy using volatility forecasts (no transaction costs)
# ----------------------------------------------------------------------

straddle_data <- var_data %>% select(date, futures_ret, sigma_model, sigma_bench)

# Realized absolute return
straddle_data <- straddle_data %>%
  mutate(abs_ret = abs(futures_ret))

# ---- FIX 1: Premium and capital from in‑sample only (no look‑ahead) ----
# In‑sample = all data before test start
test_start_date <- test_dates[1]
insample_abs_ret <- zhang_data %>%
  filter(date < test_start_date) %>%
  pull(futures_ret) %>%
  abs()

if (length(insample_abs_ret) == 0) {
  stop("No in‑sample data to compute straddle premium. Check test start date.")
}
premium <- mean(insample_abs_ret, na.rm = TRUE)
capital <- premium   # amount paid for the straddle each month

# ---- FIX 2: Use in‑sample realised volatility for thresholds (no look‑ahead) ----
# In‑sample indices: all observations before test start
insample_idx <- which(dates < test_start_date)
insample_sigma <- exp(data$LV[insample_idx] / 2)
threshold_model <- median(insample_sigma, na.rm = TRUE)
threshold_bench <- median(insample_sigma, na.rm = TRUE)

# Trading signals (1 = long straddle, -1 = short straddle)
straddle_data <- straddle_data %>%
  mutate(
    signal_model = ifelse(sigma_model > threshold_model, 1, -1),
    signal_bench = ifelse(sigma_bench > threshold_bench, 1, -1)
  )

# Payoffs (no transaction costs)
straddle_data <- straddle_data %>%
  mutate(
    payoff_model = signal_model * (abs_ret - premium),
    payoff_bench = signal_bench * (abs_ret - premium)
  )

# Cumulative profits
straddle_data <- straddle_data %>%
  mutate(
    cum_profit_model = cumsum(payoff_model),
    cum_profit_bench = cumsum(payoff_bench)
  )

# ---- FIX 3: Risk‑free rate conversion (TB3MS is annualized %) ----
zhang_rf <- read.csv("ZhangData32_3.0.csv") %>%
  select(date, TB3MS) %>%
  mutate(date = as.Date(date),
         rf_monthly = (1 + TB3MS/100)^(1/12) - 1)   # convert annual % to monthly decimal

# Merge risk‑free rate into straddle_data
straddle_data <- straddle_data %>%
  left_join(zhang_rf, by = "date") %>%
  filter(complete.cases(.))

# Compute excess returns (net of risk‑free rate)
excess_model <- (straddle_data$payoff_model / capital) - straddle_data$rf_monthly
excess_bench <- (straddle_data$payoff_bench / capital) - straddle_data$rf_monthly

# Annualised Sharpe ratios (assuming i.i.d. monthly returns)
sharpe_model <- mean(excess_model, na.rm = TRUE) / sd(excess_model, na.rm = TRUE) * sqrt(12)
sharpe_bench <- mean(excess_bench, na.rm = TRUE) / sd(excess_bench, na.rm = TRUE) * sqrt(12)

# Total profits (sum of payoffs)
total_return_model <- sum(straddle_data$payoff_model, na.rm = TRUE)
total_return_bench <- sum(straddle_data$payoff_bench, na.rm = TRUE)

cat("\n--- Straddle strategy results (no transaction costs) ---\n")
cat("Total profit (OCMT):           ", round(total_return_model, 4), "\n")
cat("Total profit (Benchmark):      ", round(total_return_bench, 4), "\n")
cat("Annualized Sharpe (OCMT):      ", round(sharpe_model, 4), "\n")
cat("Annualized Sharpe (Benchmark): ", round(sharpe_bench, 4), "\n")

# Plot cumulative profit
p_straddle <- ggplot(straddle_data, aes(x = date)) +
  geom_line(aes(y = cum_profit_model, color = "OCMT"), linewidth = 0.8) +
  geom_line(aes(y = cum_profit_bench, color = "Benchmark"), linewidth = 0.8) +
  labs(title = "Cumulative straddle profit (OCMT vs Historical mean)",
       x = "Time", y = "Cumulative profit") +
  scale_color_manual(name = "Model", values = c("OCMT" = "blue", "Benchmark" = "red")) +
  theme_classic() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_straddle)

# Optional: switching frequency
switch_model <- mean(abs(diff(straddle_data$signal_model)), na.rm = TRUE)
switch_bench <- mean(abs(diff(straddle_data$signal_bench)), na.rm = TRUE)
cat("\nAverage monthly signal change (OCMT): ", round(switch_model, 4))
cat("\nAverage monthly signal change (Benchmark): ", round(switch_bench, 4), "\n")

cat("\n=== Script finished ===\n")