

library(dplyr)
library(tidyr)
library(ggplot2)
library(rugarch)

# ----------------------------------------------------------------------
# 0. User settings
# ----------------------------------------------------------------------
DATA_FILE <- "VolatilityData.csv"
DATE_VAR  <- "date"
TARGET_VAR <- "LV"

# Expanding window 
pval_grid   <- c(0.01, 0.05, 0.10)
delta1_grid <- 1
delta2_grid <- 2          

# Split dates based on volatility sample (2008-01 to 2024-02)
TRAIN_END <- as.Date("2010-12-01")
VAL_END   <- as.Date("2013-12-01")
TEST_END  <- as.Date("2024-02-01")

# ----------------------------------------------------------------------
# 1. Helper functions
# ----------------------------------------------------------------------

# Scale training and test data using training means and standard deviations
scale_train_test <- function(X_train, X_test) {
  X_train <- as.matrix(X_train)
  X_test  <- as.matrix(X_test)
  mu <- colMeans(X_train, na.rm = TRUE)
  sdv <- apply(X_train, 2, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  X_train_sc <- sweep(X_train, 2, mu, "-")
  X_train_sc <- sweep(X_train_sc, 2, sdv, "/")
  X_test_sc <- sweep(X_test, 2, mu, "-")
  X_test_sc <- sweep(X_test_sc, 2, sdv, "/")
  list(X_train = X_train_sc, X_test = X_test_sc)
}

# Clark-West test 
clark_west_test <- function(y, f_bench, f_model) {
  d <- (y - f_bench)^2 - ((y - f_model)^2 - (f_bench - f_model)^2)
  n <- length(d)
  t_stat <- mean(d) / (sd(d) / sqrt(n))
  p_val <- 1 - pnorm(t_stat)
  list(stat = t_stat, p_value = p_val)
}

# Absolute t-statistic for a candidate variable conditional on already selected set
abs_t_first <- function(y, x, Xsel = NULL) {
  x <- as.numeric(x)
  n <- length(y)
  if (is.null(Xsel)) {
    Xsel <- matrix(nrow = n, ncol = 0)
  } else {
    Xsel <- as.matrix(Xsel)
    if (ncol(Xsel) == 0) Xsel <- matrix(nrow = n, ncol = 0)
  }
  fit <- tryCatch({
    if (ncol(Xsel) == 0) {
      df <- data.frame(y = y, x = x)
      lm(y ~ x, data = df)
    } else {
      ctrl <- as.data.frame(Xsel)
      names(ctrl) <- paste0("z", seq_len(ncol(ctrl)))
      df <- data.frame(y = y, x = x, ctrl)
      lm(y ~ ., data = df)
    }
  }, error = function(e) NULL)
  if (is.null(fit)) return(-Inf)
  sm <- summary(fit)
  rn <- rownames(sm$coefficients)
  if (!("x" %in% rn)) return(-Inf)
  abs(sm$coefficients["x", "t value"])
}

# BMT selection algorithm
boosting_glm <- function(y, X, pval, delta1, delta2) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (nrow(X) != length(y)) stop("length(y) must equal nrow(X).")
  N <- ncol(X)
  if (N < 1) stop("X must have at least one column.")
  
  p1 <- pval / (N^(delta1 - 1))
  t1 <- qnorm(1 - p1 / (2 * N))
  p2 <- pval / (N^(delta2 - 1))
  t2 <- qnorm(1 - p2 / (2 * N))
  
  ind <- rep(FALSE, N)
  
  # Initial stage: unconditional t-statistics
  ts <- rep(-Inf, N)
  for (i in seq_len(N)) {
    ts[i] <- abs_t_first(y = y, x = X[, i, drop = TRUE], Xsel = NULL)
  }
  i_max <- which.max(ts)
  if (is.finite(ts[i_max]) && ts[i_max] > t1) {
    ind[i_max] <- TRUE
  }
  
  # Boosting loop: add one variable at a time if conditional t > t2
  repeat {
    ts <- rep(-Inf, N)
    for (j in seq_len(N)) {
      if (!ind[j]) {
        ts[j] <- abs_t_first(y = y, x = X[, j, drop = TRUE],
                             Xsel = X[, ind, drop = FALSE])
      }
    }
    i_max <- which.max(ts)
    if (is.finite(ts[i_max]) && ts[i_max] > t2) {
      ind[i_max] <- TRUE
    } else {
      break
    }
  }
  ind
}

# One-step-ahead forecast using selected variables
forecast_one_step <- function(y_train, X_train, X_test, selected) {
  if (!any(selected)) return(mean(y_train, na.rm = TRUE))
  Xsel_train <- as.data.frame(X_train[, selected, drop = FALSE])
  Xsel_test  <- as.data.frame(t(X_test[selected]))
  if (ncol(Xsel_train) == 0) return(mean(y_train, na.rm = TRUE))
  names(Xsel_train) <- names(Xsel_test) <- colnames(X_train)[selected]
  fit <- tryCatch(lm(y_train ~ ., data = data.frame(y_train = y_train, Xsel_train)),
                  error = function(e) NULL)
  if (is.null(fit)) return(mean(y_train, na.rm = TRUE))
  pred <- tryCatch(as.numeric(predict(fit, newdata = Xsel_test))[1],
                   error = function(e) NA_real_)
  if (!is.finite(pred)) mean(y_train, na.rm = TRUE) else pred
}

# Recursive evaluation (expanding window)
run_forecast <- function(start_idx, end_idx, y_all, X_all, pval, d1, d2) {
  f_model <- numeric(0)
  f_bench <- numeric(0)
  y_oos   <- numeric(0)
  sel_mat <- matrix(0L, nrow = end_idx - start_idx + 1, ncol = ncol(X_all))
  colnames(sel_mat) <- colnames(X_all)
  row_i <- 1
  
  for (t in start_idx:end_idx) {
    y_train <- y_all[1:(t - 1)]
    X_train <- X_all[1:(t - 1), , drop = FALSE]
    X_test  <- X_all[t, , drop = FALSE]
    
    sc <- scale_train_test(X_train, X_test)
    X_train_sc <- sc$X_train
    X_test_sc  <- sc$X_test
    
    selected <- boosting_glm(y = y_train, X = X_train_sc,
                             pval = pval, delta1 = d1, delta2 = d2)
    
    f0 <- mean(y_train, na.rm = TRUE)      # historical mean benchmark
    f1 <- forecast_one_step(y_train, X_train_sc, X_test_sc, selected)
    
    f_model <- c(f_model, f1)
    f_bench <- c(f_bench, f0)
    y_oos   <- c(y_oos, y_all[t])
    sel_mat[row_i, ] <- as.integer(selected)
    row_i <- row_i + 1
  }
  list(
    r2 = 1 - sum((y_oos - f_model)^2) / sum((y_oos - f_bench)^2),
    f_model = f_model,
    f_bench = f_bench,
    y = y_oos,
    sel = sel_mat
  )
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
  arrange(.data[[DATE_VAR]]) %>%
  filter(complete.cases(.))

dates <- df[[DATE_VAR]]
y_all <- df[[TARGET_VAR]]
X_all <- as.matrix(df[, predictors, drop = FALSE])
colnames(X_all) <- predictors

train_end_idx <- max(which(dates <= TRAIN_END))
val_end_idx   <- max(which(dates <= VAL_END))
test_end_idx  <- max(which(dates <= TEST_END))

if (is.infinite(train_end_idx) || is.na(train_end_idx)) stop("TRAIN_END not found.")
if (is.infinite(val_end_idx) || is.na(val_end_idx)) stop("VAL_END not found.")
if (is.infinite(test_end_idx) || is.na(test_end_idx)) stop("TEST_END not found.")

cat("Training   :", as.character(min(dates[1:train_end_idx])), "to", as.character(max(dates[1:train_end_idx])), "\n")
cat("Validation :", as.character(dates[train_end_idx+1]), "to", as.character(dates[val_end_idx]), "\n")
cat("Test       :", as.character(dates[val_end_idx+1]), "to", as.character(dates[test_end_idx]), "\n")

# ----------------------------------------------------------------------
# 3. Grid search on validation period 
# ----------------------------------------------------------------------
grid <- expand.grid(pval = pval_grid, d1 = delta1_grid, d2 = delta2_grid,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$r2 <- NA_real_

for (g in seq_len(nrow(grid))) {
  res_val <- run_forecast(start_idx = train_end_idx + 1,
                          end_idx   = val_end_idx,
                          y_all     = y_all,
                          X_all     = X_all,
                          pval      = grid$pval[g],
                          d1        = grid$d1[g],
                          d2        = grid$d2[g])
  grid$r2[g] <- res_val$r2
  cat("grid", g, "of", nrow(grid), "done; R2 =", grid$r2[g], "\n")
}

best <- grid[which.max(grid$r2), ]
cat("\nBest hyperparameters:\n")
print(best)

# ----------------------------------------------------------------------
# 4. Final test evaluation
# ----------------------------------------------------------------------
final <- run_forecast(start_idx = val_end_idx + 1,
                      end_idx   = test_end_idx,
                      y_all     = y_all,
                      X_all     = X_all,
                      pval      = best$pval,
                      d1        = best$d1,
                      d2        = best$d2)

actual <- final$y
pred_bmt <- final$f_model
pred_bench <- final$f_bench

# Point forecast metrics
oos_r2 <- final$r2
mspe_bmt <- mean((actual - pred_bmt)^2)
mspe_bench <- mean((actual - pred_bench)^2)

# Directional metrics on changes
actual_diff <- diff(actual)
pred_diff <- diff(pred_bmt)
success <- mean(sign(actual_diff) == sign(pred_diff))
pt <- DACTest(pred_diff, actual_diff, test = "PT")

# Clark-West test (nested: historical mean benchmark)
cw <- clark_west_test(actual, pred_bench, pred_bmt)

cat("\n========== FINAL TEST RESULTS ==========\n")
cat("OOS R2:", oos_r2, "\n")
cat("MSPE (BMT):", mspe_bmt, "\n")
cat("MSPE (Historical mean):", mspe_bench, "\n")
cat("Clark-West stat:", cw$stat, "  p-value:", cw$p_value, "\n")
cat("Success ratio (on changes):", success, "\n")
cat("Pesaran-Timmermann test (on changes):\n")
print(pt)

# ----------------------------------------------------------------------
# 5. CSPE difference plot (BMT vs Historical mean)
# ----------------------------------------------------------------------
cspe_diff <- cumsum((actual - pred_bench)^2 - (actual - pred_bmt)^2)
test_dates <- dates[(val_end_idx + 1):test_end_idx]
cspe_df <- data.frame(Date = test_dates, CSPE_Diff = cspe_diff)

p_cspe <- ggplot(cspe_df, aes(x = Date, y = CSPE_Diff)) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  theme_classic() +
  labs(title = paste("Cumulative squared prediction error difference (BMT vs Historical mean)\n",
                     "pval =", best$pval, ", δ1 =", 1, ", δ2 =", 2),
       x = "Time", y = "CSPE (Benchmark - BMT)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cspe)

# ----------------------------------------------------------------------
# 6. Model size over time
# ----------------------------------------------------------------------
model_size <- rowSums(final$sel)
model_size_df <- data.frame(Date = test_dates, ModelSize = model_size)
p_size <- ggplot(model_size_df, aes(x = Date, y = ModelSize)) +
  geom_line(color = "darkred", linewidth = 0.8) +
  geom_point(size = 0.8, color = "darkred") +
  theme_classic() +
  labs(title = "BMT model size over time",
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
# 7. Variable importance (selection frequencies)
# ----------------------------------------------------------------------
importance <- colMeans(final$sel)
importance_sorted <- sort(importance, decreasing = TRUE)
cat("\n--- Variable selection frequencies (%) ---\n")
print(round(100 * importance_sorted, 1))

# Heatmap (ordered by frequency)
ordered_vars <- names(importance_sorted)
if (length(ordered_vars) > 0) {
  sel_plot <- final$sel[, ordered_vars, drop = FALSE]
  heat_df <- as.data.frame(sel_plot)
  heat_df$time <- test_dates
  long_df <- heat_df %>%
    pivot_longer(cols = all_of(ordered_vars),
                 names_to = "variable",
                 values_to = "selected")
  p_heat <- ggplot(long_df, aes(x = time, y = variable, fill = factor(selected))) +
    geom_tile() +
    scale_fill_manual(values = c("0" = "white", "1" = "blue"),
                      labels = c("Not selected", "Selected"), name = "") +
    labs(title = "BMT Selection Heatmap", x = "Time", y = NULL) +
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

# Align forecasts with test period dates
var_data <- data.frame(date = test_dates,
                       LV_forecast = pred_bmt,
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

cat("\n--- VaR savings (95% confidence, FHS) – BMT ---\n")
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
  labs(title = paste("VaR savings (BMT vs Historical mean) - FHS\n"),
       x = "Time", y = "VaR reduction") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_var)

# ----------------------------------------------------------------------
# Straddle strategy using volatility forecasts (no transaction costs) – BMT
# ----------------------------------------------------------------------

straddle_data <- var_data %>% select(date, futures_ret, sigma_model, sigma_bench)

# Realized absolute return
straddle_data <- straddle_data %>%
  mutate(abs_ret = abs(futures_ret))

# ---- FIX 1: Premium and capital from in‑sample only (no look‑ahead) ----
# In‑sample = all data before test start (i.e., dates < test_dates[1])
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

cat("\n--- Straddle strategy results (no transaction costs) – BMT ---\n")
cat("Total profit (BMT):       ", round(total_return_model, 4), "\n")
cat("Total profit (Benchmark): ", round(total_return_bench, 4), "\n")
cat("Annualized Sharpe (BMT):  ", round(sharpe_model, 4), "\n")
cat("Annualized Sharpe (Benchmark):", round(sharpe_bench, 4), "\n")

# Plot cumulative profit
p_straddle <- ggplot(straddle_data, aes(x = date)) +
  geom_line(aes(y = cum_profit_model, color = "BMT"), linewidth = 0.8) +
  geom_line(aes(y = cum_profit_bench, color = "Benchmark"), linewidth = 0.8) +
  labs(title = "Cumulative straddle profit (BMT vs Historical mean)",
       x = "Time", y = "Cumulative profit") +
  scale_color_manual(name = "Model", values = c("BMT" = "blue", "Benchmark" = "red")) +
  theme_classic() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_straddle)

# Optional switching frequency
switch_model <- mean(abs(diff(straddle_data$signal_model)), na.rm = TRUE)
switch_bench <- mean(abs(diff(straddle_data$signal_bench)), na.rm = TRUE)
cat("\nAverage monthly signal change (BMT): ", round(switch_model, 4))
cat("\nAverage monthly signal change (Benchmark): ", round(switch_bench, 4), "\n")


# ----------------------------------------------------------------------
# Time series plot: test period only (2014‑2024)
# ----------------------------------------------------------------------
training_end_date <- dates[val_end_idx]               # December 2013
dates_test <- dates[(val_end_idx + 1):test_end_idx]  # January 2014 – February 2024
actual_test <- final$y
bmt_test <- final$f_model
bench_test <- final$f_bench

# Full series of actual LV (light gray background)
full_actual <- data.frame(Date = dates, Actual = y_all)

plot_test_df <- data.frame(
  Date = dates_test,
  Actual = actual_test,
  BMT = bmt_test,
  Benchmark = bench_test
)

p_vol_ts <- ggplot() +
  geom_line(data = full_actual, aes(x = Date, y = Actual),
            color = "gray70", linewidth = 0.5, alpha = 0.8) +
  geom_line(data = plot_test_df, aes(x = Date, y = Actual),
            color = "black", linewidth = 0.8) +
  geom_line(data = plot_test_df, aes(x = Date, y = BMT, color = "BMT"),
            linewidth = 0.8, linetype = "dashed") +
  geom_line(data = plot_test_df, aes(x = Date, y = Benchmark, color = "Historical mean"),
            linewidth = 0.8, linetype = "dotted") +
  geom_vline(xintercept = as.numeric(training_end_date), linetype = "solid",
             color = "darkred", linewidth = 0.6) +
  annotate("text", x = training_end_date,
           y = max(y_all, na.rm = TRUE),
           label = "← Training end", hjust = -0.1, vjust = 1,
           color = "darkred", size = 3.5) +
  scale_color_manual(name = "Forecast",
                     values = c("BMT" = "blue", "Historical mean" = "forestgreen")) +
  labs(title = "Log realised volatility: actual vs. forecasts (BMT)",
       subtitle = paste("Test period (2014-2024), BMT (pval =", best$pval, ", δ1 = 1, δ2 = 2)"),
       x = "Date", y = "Log realised volatility (LV)") +
  theme_classic() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p_vol_ts)

