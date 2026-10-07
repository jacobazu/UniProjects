



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
data <- read.csv("ZhangData32_3.0.csv")
data$date <- as.Date(data$date)

# Zhang sample: Feb 1986 – Dec 2016
data <- data %>%
  filter(date >= as.Date("1986-02-01") & date <= as.Date("2025-08-01"))

# ----------------------------------------------------------------------
# 2) Define 
# ----------------------------------------------------------------------
tech_vars <- c("MA_1_9", "MA_1_12", "MA_2_9", "MA_2_12", "MA_3_9", "MA_3_12",
               "MOM_1", "MOM_2", "MOM_3", "MOM_6", "MOM_9", "MOM_12",
               "VOL_1_9", "VOL_1_12", "VOL_2_9", "VOL_2_12", "VOL_3_9", "VOL_3_12")

macro_vars <- c("TB3MS", "GS10", "infl_m", "SVOL", "epu", "kilian",
                "prod_growth", "stocks_growth", "imports_growth", "m2_growth",
                "ip_growth", "unemp_diff", "cfnai", "TCU", "BAA", "AAA")

predictors <- c(tech_vars, macro_vars)

# ----------------------------------------------------------------------
# 3) Response and other series
# ----------------------------------------------------------------------
response_col <- "target"          # spot return (fractional)
rf_col <- "TB3MS"                 # risk‑free rate (annual %)
futures_col <- "Price_Futures"    # futures price

# ----------------------------------------------------------------------
# 4) Build modelling data frame
# ----------------------------------------------------------------------
keep_cols <- c("date", response_col, rf_col, futures_col, predictors)
df <- data %>%
  select(all_of(keep_cols)) %>%
  filter(complete.cases(.))

# Compute futures simple return
df <- df %>%
  arrange(date) %>%
  mutate(futures_ret = (lead(!!sym(futures_col)) - !!sym(futures_col)) /
           !!sym(futures_col)) %>%
  filter(!is.na(futures_ret))   # drop last month 

# Extract target, predictors, dates
y <- as.numeric(df[[response_col]])
dates <- df$date
X <- as.matrix(df[, predictors])
storage.mode(X) <- "double"

cat("Sample starts at:", as.character(min(dates)), "\n")
cat("Sample ends at   :", as.character(max(dates)), "\n")
cat("N observations   :", nrow(df), "\n")

# ----------------------------------------------------------------------
# 5) Out‑of‑sample start (Jan 2001)
# ----------------------------------------------------------------------
oos_start <- which(dates >= as.Date("2008-01-01"))[1]
if (is.na(oos_start)) stop("OOS start date not found.")
cat("OOS starts at    :", as.character(dates[oos_start]), "\n")
cat("λ selection method:", lambda_method, "\n")

oos_idx <- oos_start:length(y)
in_sample_idx <- 1:(oos_start - 1)

# ----------------------------------------------------------------------
# 6) Helper functions
# ----------------------------------------------------------------------
r2_os_pct <- function(actual, model, bench) {
  mspe_model <- mean((actual - model)^2, na.rm = TRUE)
  mspe_bench <- mean((actual - bench)^2, na.rm = TRUE)
  100 * (1 - mspe_model / mspe_bench)
}

mspe <- function(actual, forecast) mean((actual - forecast)^2, na.rm = TRUE)
success_ratio_pct <- function(actual, forecast) 100 * mean(sign(actual) == sign(forecast), na.rm = TRUE)

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
# 7) In‑sample R² for ordering the heatmap
# ----------------------------------------------------------------------
insample_y <- y[in_sample_idx]
insample_r2 <- sapply(predictors, function(var) {
  fit <- lm(insample_y ~ X[in_sample_idx, var])
  summary(fit)$r.squared
})
ordered_vars_insample <- names(sort(insample_r2, decreasing = TRUE))

# ----------------------------------------------------------------------
# 8) Recursive LASSO with selected λ method
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
    cv_fit <- cv.glmnet(X_est, y_est, alpha = 1, standardize = TRUE, nfolds = 10)
    best_lambda <- cv_fit$lambda.min
  } else if (lambda_method == "BIC") {
    # Fit a path and compute BIC for each lambda
    fit_path <- glmnet(X_est, y_est, alpha = 1, standardize = TRUE)
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
# 9) Forecast evaluation (identical to before)
# ----------------------------------------------------------------------
actual_oos <- y[oos_idx]
lasso_oos <- pred_lasso[oos_idx]
bench_oos <- pred_bench[oos_idx]

eval_tbl <- data.frame(
  Statistic = c("MSPE (LASSO)", "MSPE (Historical mean)", "OOS R2 (%)", "Success ratio (%)"),
  Value = c(
    mspe(actual_oos, lasso_oos),
    mspe(actual_oos, bench_oos),
    r2_os_pct(actual_oos, lasso_oos, bench_oos),
    success_ratio_pct(actual_oos, lasso_oos)
  )
)
cat("\n--- Forecast evaluation ---\n")
print(eval_tbl, row.names = FALSE)

cw <- cw_test(actual_oos, bench_oos, lasso_oos)
cat("\nClark-West test:\n  t-stat:", round(cw["t_stat"], 4),
    "  p-value:", round(cw["p_value"], 4), "\n")

pt <- DACTest(lasso_oos, actual_oos, test = "PT")
cat("\nPesaran-Timmermann test:\n")
print(pt)

# ----------------------------------------------------------------------
# 10) CSPE difference plot
# ----------------------------------------------------------------------
cspe_diff <- cumsum((actual_oos - bench_oos)^2 - (actual_oos - lasso_oos)^2)
cspe_df <- data.frame(Date = dates[oos_idx], CSPE_Diff = cspe_diff)
p_cspe <- ggplot(cspe_df, aes(x = Date, y = CSPE_Diff)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = paste("Cumulative squared prediction error difference (", lambda_method, ")", sep = ""),
       x = "Time", y = "CSPE (Benchmark - LASSO)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cspe)

# ----------------------------------------------------------------------
# 11) Heatmap ordered by in‑sample R²
# ----------------------------------------------------------------------
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
# 12) Selection frequencies (as percentages)
# ----------------------------------------------------------------------
selected_all <- unlist(selected_vars[oos_idx], use.names = FALSE)
selected_all <- selected_all[!is.na(selected_all) & selected_all != ""]
freq <- table(selected_all)
freq_pct <- sort(100 * freq / length(oos_idx), decreasing = TRUE)

cat("\n--- Variable selection frequencies (percentage of out-of-sample periods) ---\n")
print(round(freq_pct, 1))

# Heatmap ordered by frequency
ordered_vars_freq <- names(freq_pct)
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

# ----------------------------------------------------------------------
# 13) λ evolution plot (log scale)
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

# ----------------------------------------------------------------------
# 14) Portfolio exercise 
# ----------------------------------------------------------------------
gamma <- 3
w_lower <- -1.5
w_upper <- 1.5
vol_window <- 60

asset_raw <- df$futures_ret
rf_raw <- df[[rf_col]]
rf_month <- (1 + rf_raw / 100)^(1/12) - 1
asset_excess <- asset_raw - rf_month

mu_lasso <- lasso_oos
mu_bench <- bench_oos
asset_excess_oos <- asset_excess[oos_idx]
rf_oos <- rf_month[oos_idx]

sigma2 <- rep(NA_real_, length(actual_oos))
for (t in seq_along(oos_idx)) {
  global_t <- oos_idx[t]
  if (global_t <= vol_window) next
  hist_idx <- (global_t - vol_window):(global_t - 1)
  sigma2[t] <- var(asset_excess[hist_idx], na.rm = TRUE)
}
first_sigma <- which(!is.na(sigma2))[1]
if (!is.na(first_sigma) && first_sigma > 1) sigma2[1:(first_sigma-1)] <- sigma2[first_sigma]

valid <- is.finite(mu_lasso) & is.finite(mu_bench) & is.finite(asset_excess_oos) &
  is.finite(rf_oos) & is.finite(sigma2)
mu_lasso_p <- mu_lasso[valid]
mu_bench_p <- mu_bench[valid]
asset_excess_p <- asset_excess_oos[valid]
rf_oos_p <- rf_oos[valid]
sigma2_p <- sigma2[valid]
dates_p <- dates[oos_idx][valid]

portfolio_metrics <- function(mu_hat, asset_excess_oos, rf_oos, sigma2_vec,
                              gamma = 3, w_lower = -1.5, w_upper = 1.5) {
  w <- (1 / gamma) * (mu_hat / sigma2_vec)
  w <- pmax(pmin(w, w_upper), w_lower)
  port_excess <- w * asset_excess_oos
  port_total <- rf_oos + port_excess
  wealth <- cumprod(1 + port_total)
  cer_monthly <- mean(port_total, na.rm = TRUE) - 0.5 * gamma * var(port_total, na.rm = TRUE)
  sharpe_annual <- sqrt(12) * mean(port_excess, na.rm = TRUE) / sd(port_excess, na.rm = TRUE)
  list(port_total = port_total, wealth = wealth, cer_monthly = cer_monthly,
       sharpe_annual = sharpe_annual)
}

port_lasso <- portfolio_metrics(mu_lasso_p, asset_excess_p, rf_oos_p, sigma2_p)
port_bench <- portfolio_metrics(mu_bench_p, asset_excess_p, rf_oos_p, sigma2_p)
cer_gain_annual_pct <- 1200 * (port_lasso$cer_monthly - port_bench$cer_monthly)

cat("\n--- Portfolio exercise (futures returns) ---\n")
cat("Annualized CER gain (LASSO - Historical mean) [%] :", round(cer_gain_annual_pct, 4), "\n")
cat("Annualized Sharpe ratio (LASSO)                   :", round(port_lasso$sharpe_annual, 4), "\n")
cat("Annualized Sharpe ratio (Historical mean)         :", round(port_bench$sharpe_annual, 4), "\n")

wealth_df <- data.frame(Date = dates_p,
                        LASSO = port_lasso$wealth,
                        HistoricalMean = port_bench$wealth)
wealth_long <- melt(wealth_df, id.vars = "Date", variable.name = "Strategy", value.name = "Wealth")
p_wealth <- ggplot(wealth_long, aes(Date, Wealth, color = Strategy)) +
  geom_line(linewidth = 0.8) +
  theme_classic() +
  labs(title = paste("Cumulative wealth: LASSO vs Historical Mean (", lambda_method, ")", sep = ""),
       x = "Time", y = "Wealth (initial = 1)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_wealth)


# ----------------------------------------------------------------------
# Time series plot: actual returns, LASSO forecasts, benchmark forecasts
# ----------------------------------------------------------------------

# Training end date (last month used for estimation before first OOS forecast)
training_end_date <- dates[oos_start - 1]

# Build data frame for the out‑of‑sample period
plot_oos_df <- data.frame(
  Date = dates[oos_idx],
  Actual = y[oos_idx],
  LASSO = pred_lasso[oos_idx],
  Benchmark = pred_bench[oos_idx]
)

# Full series of actual returns (for background)
full_actual <- data.frame(Date = dates, Actual = y)

# Create the plot
p_returns_ts <- ggplot() +
  # Full actual returns (light gray, thin)
  geom_line(data = full_actual, aes(x = Date, y = Actual),
            color = "gray70", linewidth = 0.5, alpha = 0.8) +
  # Out‑of‑sample actual returns (black, thicker)
  geom_line(data = plot_oos_df, aes(x = Date, y = Actual),
            color = "black", linewidth = 0.8) +
  # LASSO forecasts
  geom_line(data = plot_oos_df, aes(x = Date, y = LASSO, color = "LASSO"),
            linewidth = 0.8, linetype = "dashed") +
  # Historical mean forecasts
  geom_line(data = plot_oos_df, aes(x = Date, y = Benchmark, color = "Historical mean"),
            linewidth = 0.8, linetype = "dotted") +
  # Vertical line at end of training period
  geom_vline(xintercept = as.numeric(training_end_date), linetype = "solid",
             color = "darkred", linewidth = 0.6) +
  # Annotation for training end
  annotate("text", x = training_end_date,
           y = max(y, na.rm = TRUE),
           label = "← Training end", hjust = -0.1, vjust = 1,
           color = "darkred", size = 3.5) +
  scale_color_manual(name = "Forecast",
                     values = c("LASSO" = "blue", "Historical mean" = "forestgreen")) +
  labs(title = paste("Oil spot returns: actual vs. forecasts (", lambda_method, ")", sep = ""),
       subtitle = paste("LASSO (", lambda_method, " λ selection), Historical mean benchmark"),
       x = "Date", y = "Spot return") +
  theme_classic() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

print(p_returns_ts)
