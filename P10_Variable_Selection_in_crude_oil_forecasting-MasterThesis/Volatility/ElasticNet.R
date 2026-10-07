

library(glmnet)
library(dplyr)
library(reshape2)
library(ggplot2)
library(forecast)  
library(rugarch)

# ----------------------------------------------------------------------
# 1) Load and prepare data
# ----------------------------------------------------------------------
data <- read.csv("VolatilityData.csv")
data$date <- as.Date(data$date)

# Ensure sorted by date
data <- data %>% arrange(date)

# ----------------------------------------------------------------------
# 2) Define predictors (all columns except date and LV)
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
# 5) Recursive Elastic Net with 60/40 split for α and λ
# ----------------------------------------------------------------------
# Grid for α (mixing parameter)
alpha_grid <- c(0.1, 0.3, 0.5, 0.7, 0.9)

n <- length(y)
pred_enet <- rep(NA_real_, n)
pred_bench <- rep(NA_real_, n)
lambda_store <- rep(NA_real_, n)
alpha_store <- rep(NA_real_, n)
selected_vars <- vector("list", n)

for (i in oos_idx) {
  est_idx <- 1:(i - 1)
  n_est <- length(est_idx)
  
  # Validation block = last 40% of estimation window
  n_val <- max(5, floor(0.40 * n_est))
  n_tr <- n_est - n_val
  if (n_tr < 10) next
  
  tr_idx <- est_idx[1:n_tr]
  val_idx <- est_idx[(n_tr + 1):n_est]
  
  X_tr <- X[tr_idx, , drop = FALSE]
  y_tr <- y[tr_idx]
  X_val <- X[val_idx, , drop = FALSE]
  y_val <- y[val_idx]
  
  # Grid search over α
  best_alpha <- NA
  best_lambda <- NA
  best_val_mse <- Inf
  
  for (alpha_val in alpha_grid) {
    # Fit elastic net path on training block
    fit_path <- tryCatch(
      glmnet(X_tr, y_tr, alpha = alpha_val, standardize = TRUE),
      error = function(e) NULL
    )
    if (is.null(fit_path)) next
    
    lambdas <- fit_path$lambda
    
    # Predict on validation set for all lambdas
    val_pred <- predict(fit_path, newx = X_val, s = lambdas)
    y_val_mat <- matrix(y_val, nrow = length(y_val), ncol = length(lambdas))
    val_mse <- colMeans((y_val_mat - val_pred)^2)
    
    # Find best lambda for this α
    best_lambda_alpha <- lambdas[which.min(val_mse)]
    best_mse_alpha <- min(val_mse)
    
    if (best_mse_alpha < best_val_mse) {
      best_val_mse <- best_mse_alpha
      best_alpha <- alpha_val
      best_lambda <- best_lambda_alpha
    }
  }
  
  if (is.na(best_alpha)) {
    # Fallback: use historical mean for this period
    pred_enet[i] <- mean(y[est_idx], na.rm = TRUE)
    pred_bench[i] <- mean(y[est_idx], na.rm = TRUE)
    lambda_store[i] <- NA
    alpha_store[i] <- NA
    selected_vars[[i]] <- character(0)
    next
  }
  
  lambda_store[i] <- best_lambda
  alpha_store[i] <- best_alpha
  
  # Refit on all data up to i-1 with best (α, λ)
  final_fit <- glmnet(X[est_idx, , drop = FALSE], y[est_idx],
                      alpha = best_alpha, lambda = best_lambda, standardize = TRUE)
  pred_enet[i] <- as.numeric(predict(final_fit, newx = X[i, , drop = FALSE]))
  pred_bench[i] <- mean(y[est_idx], na.rm = TRUE)
  
  # Store selected predictors
  beta <- as.matrix(coef(final_fit, s = best_lambda))
  selected_vars[[i]] <- predictors[which(beta[-1, 1] != 0)]
}

# ----------------------------------------------------------------------
# 6) Forecast evaluation
# ----------------------------------------------------------------------
actual_oos <- y[oos_idx]
enet_oos <- pred_enet[oos_idx]
bench_oos <- pred_bench[oos_idx]

eval_tbl <- data.frame(
  Statistic = c("MSPE (ENet)", "MSPE (Historical mean)", "OOS R2 (%)", "Success ratio (%)"),
  Value = c(
    mspe(actual_oos, enet_oos),
    mspe(actual_oos, bench_oos),
    100 * (1 - mspe(actual_oos, enet_oos) / mspe(actual_oos, bench_oos)),
    success_ratio_pct(actual_oos, enet_oos)
  )
)
cat("\n--- Forecast evaluation ---\n")
print(eval_tbl, row.names = FALSE)

cw <- cw_test(actual_oos, bench_oos, enet_oos)
cat("\nClark-West test:\n  t-stat:", round(cw["t_stat"], 4),
    "  p-value:", round(cw["p_value"], 4), "\n")

# Directional accuracy test (Pesaran–Timmermann) on changes
actual_diff <- diff(actual_oos)
forecast_diff <- diff(enet_oos)
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
cspe_diff <- cumsum((actual_oos - bench_oos)^2 - (actual_oos - enet_oos)^2)
cspe_df <- data.frame(Date = dates[oos_idx], CSPE_Diff = cspe_diff)
p_cspe <- ggplot(cspe_df, aes(x = Date, y = CSPE_Diff)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "Cumulative squared prediction error difference (ENet vs Historical mean)",
       x = "Time", y = "CSPE (Benchmark - ENet)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cspe)

# ----------------------------------------------------------------------
# 8) Selection matrix and heatmap
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
  labs(title = "ENet Predictor Selection ",
       x = "Time", y = "Predictor") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
print(p_heat_insample)

# ----------------------------------------------------------------------
# 9) Selection frequencies 
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
    labs(title = "ENet Predictor Selection (ordered by frequency)",
         x = "Time", y = "Predictor") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank())
  print(p_heat_freq)
}

# ----------------------------------------------------------------------
# 10) λ and α evolution plots
# ----------------------------------------------------------------------
lambda_oos <- lambda_store[oos_idx]
alpha_oos <- alpha_store[oos_idx]

lambda_df <- data.frame(Date = dates[oos_idx], Lambda = lambda_oos)
p_lambda <- ggplot(lambda_df, aes(x = Date, y = log(Lambda))) +
  geom_line(color = "blue", linewidth = 0.8) +
  theme_classic() +
  labs(title = "Evolution of chosen λ (log scale) – ENet",
       x = "Time", y = "log(λ)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_lambda)

alpha_df <- data.frame(Date = dates[oos_idx], Alpha = alpha_oos)
p_alpha <- ggplot(alpha_df, aes(x = Date, y = Alpha)) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  theme_classic() +
  labs(title = "Evolution of chosen α – ENet",
       x = "Time", y = "α") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_alpha)

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
                       LV_forecast = enet_oos,
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

cat("\n--- VaR savings (95% confidence, FHS) – Elastic Net ---\n")
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
  labs(title = "VaR savings (ENet vs Historical mean) - Filtered Historical Simulation",
       x = "Time", y = "VaR reduction") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_var)


# ----------------------------------------------------------------------
# Straddle strategy using volatility forecasts
# ----------------------------------------------------------------------

straddle_data <- var_data %>% select(date, futures_ret, sigma_model, sigma_bench)

# Realized absolute return
straddle_data <- straddle_data %>%
  mutate(abs_ret = abs(futures_ret))


insample_abs_ret <- zhang_data %>%
  filter(date < oos_start_date) %>%
  pull(futures_ret) %>%
  abs()

if (length(insample_abs_ret) == 0) {
  stop("No in‑sample data to compute straddle premium. Check oos_start_date.")
}
premium <- mean(insample_abs_ret, na.rm = TRUE)
capital <- premium   # amount paid for the straddle each month


insample_sigma <- exp(data$LV[in_sample_idx] / 2)
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


zhang_rf <- read.csv("ZhangData32_3.0.csv") %>%
  select(date, TB3MS) %>%
  mutate(date = as.Date(date),
         rf_monthly = (1 + TB3MS/100)^(1/12) - 1)   # convert annual % to monthly decimal

# Merge risk‑free rate into straddle_data
straddle_data <- straddle_data %>%
  left_join(zhang_rf, by = "date") %>%
  filter(complete.cases(.))

# Compute excess returns 
excess_model <- (straddle_data$payoff_model / capital) - straddle_data$rf_monthly
excess_bench <- (straddle_data$payoff_bench / capital) - straddle_data$rf_monthly

# Annualised Sharpe ratios 
sharpe_model <- mean(excess_model, na.rm = TRUE) / sd(excess_model, na.rm = TRUE) * sqrt(12)
sharpe_bench <- mean(excess_bench, na.rm = TRUE) / sd(excess_bench, na.rm = TRUE) * sqrt(12)

# Total profits (sum of payoffs)
total_return_model <- sum(straddle_data$payoff_model, na.rm = TRUE)
total_return_bench <- sum(straddle_data$payoff_bench, na.rm = TRUE)

cat("\n--- Straddle strategy results (no transaction costs) ---\n")
cat("Total profit (ENet):           ", round(total_return_model, 4), "\n")
cat("Total profit (Benchmark):      ", round(total_return_bench, 4), "\n")
cat("Annualized Sharpe (ENet):      ", round(sharpe_model, 4), "\n")
cat("Annualized Sharpe (Benchmark): ", round(sharpe_bench, 4), "\n")

# Plot cumulative profit
p_straddle <- ggplot(straddle_data, aes(x = date)) +
  geom_line(aes(y = cum_profit_model, color = "ENet"), linewidth = 0.8) +
  geom_line(aes(y = cum_profit_bench, color = "Benchmark"), linewidth = 0.8) +
  labs(title = "Cumulative straddle profit (ENet vs Historical mean)",
       x = "Time", y = "Cumulative profit") +
  scale_color_manual(name = "Model", values = c("ENet" = "blue", "Benchmark" = "red")) +
  theme_classic() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_straddle)

#switching frequency
switch_model <- mean(abs(diff(straddle_data$signal_model)), na.rm = TRUE)
switch_bench <- mean(abs(diff(straddle_data$signal_bench)), na.rm = TRUE)
cat("\nAverage monthly signal change (ENet): ", round(switch_model, 4))
cat("\nAverage monthly signal change (Benchmark): ", round(switch_bench, 4), "\n")





