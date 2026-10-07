

library(glmnet)
library(dplyr)
library(reshape2)
library(ggplot2)
library(forecast)  
library(rugarch)   

# ----------------------------------------------------------------------
# 0) Choose λ selection method: "CV" or "BIC"
# ----------------------------------------------------------------------
lambda_method <- "BIC"   

# ----------------------------------------------------------------------
# 1) Load and prepare data
# ----------------------------------------------------------------------
data <- read.csv("VolatilityData.csv")
data$date <- as.Date(data$date)

# Ensure sorted by date
data <- data %>% arrange(date)

# ----------------------------------------------------------------------
# 2) Define predictors 
# ----------------------------------------------------------------------
predictors <- setdiff(names(data), c("date", "LV"))
X <- as.matrix(data[, predictors])
y <- data$LV   # target: log realized volatility

# Remove incomplete rows
complete_idx <- complete.cases(X, y)
X <- X[complete_idx, ]
y <- y[complete_idx]
data <- data[complete_idx, ]
dates <- data$date

cat("Sample starts at:", as.character(min(dates)), "\n")
cat("Sample ends at   :", as.character(max(dates)), "\n")
cat("N observations   :", nrow(data), "\n")

# ----------------------------------------------------------------------
# 3) Out‑of‑sample start 
# ----------------------------------------------------------------------
oos_start_date <- as.Date("2014-01-01")
oos_start <- which(dates >= oos_start_date)[1]
if (is.na(oos_start)) stop("OOS start date not found.")
cat("OOS starts at    :", as.character(dates[oos_start]), "\n")
cat("λ selection method:", lambda_method, "\n")

oos_idx <- oos_start:length(y)
in_sample_idx <- 1:(oos_start - 1)

# ----------------------------------------------------------------------
# 4) Helper functions
# ----------------------------------------------------------------------
mspe <- function(actual, forecast) mean((actual - forecast)^2, na.rm = TRUE)

success_ratio_pct <- function(actual, forecast) {
  actual_diff <- diff(actual)
  forecast_diff <- diff(forecast)
  if (length(forecast_diff) == length(actual_diff)) {
    correct <- sign(actual_diff) == sign(forecast_diff)
    return(100 * mean(correct, na.rm = TRUE))
  } else {
    return(100 * mean(sign(actual) == sign(forecast), na.rm = TRUE))
  }
}

cw_test <- function(actual, bench, model) {
  e0 <- actual - bench
  e1 <- actual - model
  d <- e0^2 - (e1^2 - (bench - model)^2)
  fit <- lm(d ~ 1)
  t_stat <- coef(summary(fit))[1, 3]
  p_val <- pt(t_stat, df = length(d) - 1, lower.tail = FALSE)
  c(t_stat = unname(t_stat), p_value = unname(p_val))
}

# ----------------------------------------------------------------------
# 5) Recursive LASSO with CV or BIC for λ
# ----------------------------------------------------------------------
n <- length(y)
pred_lasso <- rep(NA_real_, n)
pred_bench <- rep(NA_real_, n)
lambda_store <- rep(NA_real_, n)
selected_vars <- vector("list", n)

for (i in oos_idx) {
  est_idx <- 1:(i - 1)
  X_est <- X[est_idx, , drop = FALSE]
  y_est <- y[est_idx]
  
  if (lambda_method == "CV") {
    # 10‑fold cross‑validation
    cv_fit <- tryCatch(
      cv.glmnet(X_est, y_est, alpha = 1, standardize = TRUE, nfolds = 10),
      error = function(e) NULL
    )
    if (is.null(cv_fit)) next
    best_lambda <- cv_fit$lambda.min
  } else if (lambda_method == "BIC") {
    # Fit a path and compute BIC for each lambda
    fit_path <- tryCatch(
      glmnet(X_est, y_est, alpha = 1, standardize = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit_path)) next
    lambdas <- fit_path$lambda
    bic_vals <- sapply(lambdas, function(lam) {
      fit <- glmnet(X_est, y_est, alpha = 1, lambda = lam, standardize = TRUE)
      coefs <- as.matrix(coef(fit))
      n_obs <- length(y_est)
      rss <- sum((y_est - predict(fit, newx = X_est))^2)
      df <- sum(coefs[-1] != 0)
      n_obs * log(rss / n_obs) + df * log(n_obs)
    })
    best_lambda <- lambdas[which.min(bic_vals)]
  } else {
    stop("lambda_method must be 'CV' or 'BIC'")
  }
  
  lambda_store[i] <- best_lambda
  
  # Refit on all estimation data with chosen λ
  final_fit <- glmnet(X_est, y_est, alpha = 1, lambda = best_lambda, standardize = TRUE)
  pred_lasso[i] <- as.numeric(predict(final_fit, newx = X[i, , drop = FALSE]))
  pred_bench[i] <- mean(y_est, na.rm = TRUE)
  
  # Store selected predictors
  beta <- as.matrix(coef(final_fit, s = best_lambda))
  selected_vars[[i]] <- predictors[which(beta[-1, 1] != 0)]
}

# ----------------------------------------------------------------------
# 6) Forecast evaluation
# ----------------------------------------------------------------------
actual_oos <- y[oos_idx]
lasso_oos <- pred_lasso[oos_idx]
bench_oos <- pred_bench[oos_idx]

eval_tbl <- data.frame(
  Statistic = c("MSPE (LASSO)", "MSPE (Historical mean)", "OOS R2 (%)", "Success ratio (%)"),
  Value = c(
    mspe(actual_oos, lasso_oos),
    mspe(actual_oos, bench_oos),
    100 * (1 - mspe(actual_oos, lasso_oos) / mspe(actual_oos, bench_oos)),
    success_ratio_pct(actual_oos, lasso_oos)
  )
)
cat("\n--- Forecast evaluation ---\n")
print(eval_tbl, row.names = FALSE)

cw <- cw_test(actual_oos, bench_oos, lasso_oos)
cat("\nClark-West test:\n  t-stat:", round(cw["t_stat"], 4),
    "  p-value:", round(cw["p_value"], 4), "\n")

# Directional accuracy test (Pesaran–Timmermann) on changes
actual_diff <- diff(actual_oos)
forecast_diff <- diff(lasso_oos)
if (length(actual_diff) == length(forecast_diff)) {
  pt <- DACTest(forecast_diff, actual_diff, test = "PT")
  cat("\nPesaran-Timmermann test (on changes):\n")
  print(pt)
} else {
  cat("\nSkipping PT test due to length mismatch.\n")
}

# ----------------------------------------------------------------------
# 7) CSPE difference plot
# ----------------------------------------------------------------------
cspe_diff <- cumsum((actual_oos - bench_oos)^2 - (actual_oos - lasso_oos)^2)
cspe_df <- data.frame(Date = dates[oos_idx], CSPE_Diff = cspe_diff)
p_cspe <- ggplot(cspe_df, aes(x = Date, y = CSPE_Diff)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = paste("Cumulative squared prediction error difference (LASSO vs Historical mean)\nλ selection:", lambda_method),
       x = "Time", y = "CSPE (Benchmark - LASSO)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cspe)

# ----------------------------------------------------------------------
# 8) Selection matrix and heatmap (ordered by in‑sample R²)
# ----------------------------------------------------------------------
insample_y <- y[in_sample_idx]
insample_r2 <- sapply(predictors, function(var) {
  fit <- lm(insample_y ~ X[in_sample_idx, var])
  summary(fit)$r.squared
})
ordered_vars_insample <- names(sort(insample_r2, decreasing = TRUE))

selection_matrix <- matrix(0, nrow = length(predictors), ncol = length(oos_idx),
                           dimnames = list(predictors, as.character(dates[oos_idx])))
for (k in seq_along(oos_idx)) {
  ii <- oos_idx[k]
  vars_k <- selected_vars[[ii]]
  if (!is.null(vars_k) && length(vars_k) > 0) selection_matrix[vars_k, k] <- 1
}

selection_plot_insample <- selection_matrix[ordered_vars_insample, , drop = FALSE]
heat_data_insample <- melt(selection_plot_insample)
colnames(heat_data_insample) <- c("Variable", "Time", "Selected")
heat_data_insample$Time <- as.Date(heat_data_insample$Time)
heat_data_insample$Variable <- factor(heat_data_insample$Variable, levels = ordered_vars_insample)

p_heat_insample <- ggplot(heat_data_insample, aes(Time, Variable, fill = factor(Selected))) +
  geom_raster() +
  scale_fill_manual(values = c("0" = "white", "1" = "blue"),
                    name = "", labels = c("Not Selected", "Selected")) +
  scale_x_date(date_breaks = "12 months", date_labels = "%Y") +
  labs(title = paste("LASSO Predictor Selection (", lambda_method, ")", sep = ""),
       x = "Time", y = "Predictor") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
print(p_heat_insample)

# ----------------------------------------------------------------------
# 9) Selection frequencies (as percentages)
# ----------------------------------------------------------------------
selected_all <- unlist(selected_vars[oos_idx], use.names = FALSE)
selected_all <- selected_all[!is.na(selected_all) & selected_all != ""]
freq <- table(selected_all)
freq_pct <- sort(100 * freq / length(oos_idx), decreasing = TRUE)
cat("\n--- Variable selection frequencies (percentage of out-of-sample periods) ---\n")
print(round(freq_pct, 1))

# Heatmap ordered by frequency
ordered_vars_freq <- names(freq_pct)
if (length(ordered_vars_freq) > 0) {
  selection_plot_freq <- selection_matrix[ordered_vars_freq, , drop = FALSE]
  heat_data_freq <- melt(selection_plot_freq)
  colnames(heat_data_freq) <- c("Variable", "Time", "Selected")
  heat_data_freq$Time <- as.Date(heat_data_freq$Time)
  heat_data_freq$Variable <- factor(heat_data_freq$Variable, levels = ordered_vars_freq)
  
  p_heat_freq <- ggplot(heat_data_freq, aes(Time, Variable, fill = factor(Selected))) +
    geom_tile() +
    scale_fill_manual(values = c("white", "blue"), name = "",
                      labels = c("Not Selected", "Selected")) +
    scale_x_date(date_breaks = "12 months", date_labels = "%Y") +
    labs(title = paste("LASSO Predictor Selection (", lambda_method, ", ordered by frequency)", sep = ""),
         x = "Time", y = "Predictor") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  print(p_heat_freq)
}

# ----------------------------------------------------------------------
# 10) λ evolution plot (log scale)
# ----------------------------------------------------------------------
lambda_oos <- lambda_store[oos_idx]
lambda_df <- data.frame(Date = dates[oos_idx], Lambda = lambda_oos)
p_lambda <- ggplot(lambda_df, aes(x = Date, y = log(Lambda))) +
  geom_line(color = "blue", linewidth = 0.8) +
  theme_classic() +
  labs(title = paste("Evolution of chosen λ (log scale) – ", lambda_method, sep = ""),
       x = "Time", y = "log(λ)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_lambda)

cat("\n=== Script finished ===\n")



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

# Align forecasts with dates
forecast_dates <- dates[oos_idx]
var_data <- data.frame(date = forecast_dates,
                       LV_forecast = lasso_oos,
                       LV_bench = bench_oos) %>%
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
    empirical_quantile[i] <- qnorm(0.05)
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
  labs(title = paste("VaR savings (LASSO vs Historical mean) - FHS\nλ selection:", lambda_method),
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

# Thresholds: median of forecasted volatility
threshold_model <- median(straddle_data$sigma_model, na.rm = TRUE)
threshold_bench <- median(straddle_data$sigma_bench, na.rm = TRUE)

# Straddle premium: long‑term average of absolute returns
premium <- mean(straddle_data$abs_ret, na.rm = TRUE)

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

# Performance metrics
sharpe_model <- mean(straddle_data$payoff_model, na.rm = TRUE) / 
  sd(straddle_data$payoff_model, na.rm = TRUE) * sqrt(12)
sharpe_bench <- mean(straddle_data$payoff_bench, na.rm = TRUE) / 
  sd(straddle_data$payoff_bench, na.rm = TRUE) * sqrt(12)
total_return_model <- sum(straddle_data$payoff_model, na.rm = TRUE)
total_return_bench <- sum(straddle_data$payoff_bench, na.rm = TRUE)

cat("\n--- Straddle strategy results (no transaction costs) ---\n")
cat("Total profit (LASSO):       ", round(total_return_model, 4), "\n")
cat("Total profit (Benchmark):   ", round(total_return_bench, 4), "\n")
cat("Annualized Sharpe (LASSO):  ", round(sharpe_model, 4), "\n")
cat("Annualized Sharpe (Benchmark):", round(sharpe_bench, 4), "\n")

# Plot cumulative profit
p_straddle <- ggplot(straddle_data, aes(x = date)) +
  geom_line(aes(y = cum_profit_model, color = "LASSO"), linewidth = 0.8) +
  geom_line(aes(y = cum_profit_bench, color = "Benchmark"), linewidth = 0.8) +
  labs(title = paste("Cumulative straddle profit (LASSO vs Historical mean\nλ selection:", lambda_method),
       x = "Time", y = "Cumulative profit") +
  scale_color_manual(name = "Model", values = c("LASSO" = "blue", "Benchmark" = "red")) +
  theme_classic() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_straddle)

# Optional: switching frequency
switch_model <- mean(abs(diff(straddle_data$signal_model)), na.rm = TRUE)
switch_bench <- mean(abs(diff(straddle_data$signal_bench)), na.rm = TRUE)
cat("\nAverage monthly signal change (LASSO): ", round(switch_model, 4))
cat("\nAverage monthly signal change (Benchmark): ", round(switch_bench, 4), "\n")





# ----------------------------------------------------------------------
# Time series plot: actual LV, LASSO forecasts, historical mean benchmark
# ----------------------------------------------------------------------

# Training end date (last month before first OOS forecast)
training_end_date <- dates[oos_start - 1]   # December 2013

# Out-of-sample period dates and values
oos_dates <- dates[oos_idx]
actual_oos_plot <- y[oos_idx]
lasso_oos_plot <- pred_lasso[oos_idx]
bench_oos_plot <- pred_bench[oos_idx]

# Full series of actual LV (for background)
full_actual <- data.frame(Date = dates, Actual = y)

# Data frame for OOS period
plot_oos_df <- data.frame(
  Date = oos_dates,
  Actual = actual_oos_plot,
  LASSO = lasso_oos_plot,
  Benchmark = bench_oos_plot
)

# Create the plot
p_vol_ts <- ggplot() +
  # Full actual LV (light gray, thin)
  geom_line(data = full_actual, aes(x = Date, y = Actual),
            color = "gray70", linewidth = 0.5, alpha = 0.8) +
  # Out-of-sample actual LV (black, thicker)
  geom_line(data = plot_oos_df, aes(x = Date, y = Actual),
            color = "black", linewidth = 0.8) +
  # LASSO forecasts
  geom_line(data = plot_oos_df, aes(x = Date, y = LASSO, color = "LASSO"),
            linewidth = 0.8, linetype = "dashed") +
  # Historical mean forecasts
  geom_line(data = plot_oos_df, aes(x = Date, y = Benchmark, color = "Historical mean"),
            linewidth = 0.8, linetype = "dotted") +
  # Vertical line at end of training period (Dec 2013)
  geom_vline(xintercept = as.numeric(training_end_date), linetype = "solid",
             color = "darkred", linewidth = 0.6) +
  # Annotation for training end
  annotate("text", x = training_end_date,
           y = max(y, na.rm = TRUE),
           label = "← Training end", hjust = -0.1, vjust = 1,
           color = "darkred", size = 3.5) +
  scale_color_manual(name = "Forecast",
                     values = c("LASSO" = "blue", "Historical mean" = "forestgreen")) +
  labs(title = paste("Log realised volatility: actual vs. forecasts (LASSO, λ =", lambda_method, ")"),
       subtitle = "LASSO with expanding window, historical mean benchmark",
       x = "Date", y = "Log realised volatility (LV)") +
  theme_classic() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p_vol_ts)

