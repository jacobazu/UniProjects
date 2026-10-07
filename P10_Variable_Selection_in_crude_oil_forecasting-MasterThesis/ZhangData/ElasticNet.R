

library(glmnet)
library(dplyr)
library(reshape2)
library(ggplot2)
library(forecast)  
library(rugarch)

# ----------------------------------------------------------------------
# 1) Load and prepare data
# ----------------------------------------------------------------------
data <- read.csv("ZhangData32_3.0.csv")
data$date <- as.Date(data$date)

# Zhang sample: Feb 1986 – Dec 2016
data <- data %>%
  filter(date >= as.Date("1986-02-01") & date <= as.Date("2025-08-01"))

# ----------------------------------------------------------------------
# 2) Define predictors (all available in the CSV)
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
response_col <- "target"  
rf_col <- "TB3MS"             
futures_col <- "Price_Futures"   

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
# 8) Recursive Elastic Net with 60/40 split for α and λ
# ----------------------------------------------------------------------

alpha_grid <- seq(0, 0.9, by = 0.1) 

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
  
  # We will evaluate each alpha and choose the best (alpha, lambda) pair
  best_mse <- Inf
  best_alpha <- NA
  best_lambda <- NA
  best_fit <- NULL  
  
  for (a in alpha_grid) {
    # Fit path on training block with current alpha
    fit_path <- glmnet(X_tr, y_tr, alpha = a, standardize = TRUE)
    lambdas <- fit_path$lambda
    
    # Validate each lambda
    val_pred <- predict(fit_path, newx = X_val, s = lambdas)
    y_val_mat <- matrix(y_val, nrow = length(y_val), ncol = length(lambdas))
    val_mse <- colMeans((y_val_mat - val_pred)^2)
    
    # Best lambda for this alpha
    best_lambda_this <- lambdas[which.min(val_mse)]
    min_mse_this <- min(val_mse)
    
    if (min_mse_this < best_mse) {
      best_mse <- min_mse_this
      best_alpha <- a
      best_lambda <- best_lambda_this
    }
  }
  
  # Store selected alpha and lambda
  alpha_store[i] <- best_alpha
  lambda_store[i] <- best_lambda
  
  # Refit on all estimation data with the best (alpha, lambda)
  final_fit <- glmnet(X[est_idx, , drop = FALSE], y[est_idx],
                      alpha = best_alpha, lambda = best_lambda, standardize = TRUE)
  pred_enet[i] <- as.numeric(predict(final_fit, newx = X[i, , drop = FALSE]))
  pred_bench[i] <- mean(y[est_idx], na.rm = TRUE)
  
  # Store selected predictors
  beta <- as.matrix(coef(final_fit, s = best_lambda))
  selected_vars[[i]] <- predictors[which(beta[-1, 1] != 0)]
}

# ----------------------------------------------------------------------
# 9) Forecast evaluation
# ----------------------------------------------------------------------
actual_oos <- y[oos_idx]
enet_oos <- pred_enet[oos_idx]
bench_oos <- pred_bench[oos_idx]

eval_tbl <- data.frame(
  Statistic = c("MSPE (ENet)", "MSPE (Historical mean)", "OOS R2 (%)", "Success ratio (%)"),
  Value = c(
    mspe(actual_oos, enet_oos),
    mspe(actual_oos, bench_oos),
    r2_os_pct(actual_oos, enet_oos, bench_oos),
    success_ratio_pct(actual_oos, enet_oos)
  )
)
cat("\n--- Forecast evaluation ---\n")
print(eval_tbl, row.names = FALSE)

cw <- cw_test(actual_oos, bench_oos, enet_oos)
cat("\nClark-West test:\n  t-stat:", round(cw["t_stat"], 4),
    "  p-value:", round(cw["p_value"], 4), "\n")

# Directional accuracy test (Pesaran–Timmermann)
pt <- DACTest(enet_oos, actual_oos, test = "PT")
cat("\nPesaran-Timmermann test:\n")
print(pt)

# ----------------------------------------------------------------------
# 10) CSPE difference plot
# ----------------------------------------------------------------------
cspe_diff <- cumsum((actual_oos - bench_oos)^2 - (actual_oos - enet_oos)^2)
cspe_df <- data.frame(Date = dates[oos_idx], CSPE_Diff = cspe_diff)
p_cspe <- ggplot(cspe_df, aes(x = Date, y = CSPE_Diff)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(title = "Cumulative squared prediction error difference (ENet)",
       x = "Time", y = "CSPE (Benchmark - ENet)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cspe)

# ----------------------------------------------------------------------
# 11) Heatmap ordered by in‑sample R² (Zhang style)
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
  geom_tile() +
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
  labs(title = "ENet Predictor Selection (ordered by frequency)",
       x = "Time", y = "Predictor") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())
print(p_heat_freq)

# ----------------------------------------------------------------------
# 13) λ and α evolution plots (log scale)
# ----------------------------------------------------------------------
lambda_oos <- lambda_store[oos_idx]
alpha_oos <- alpha_store[oos_idx]

# λ evolution
lambda_df <- data.frame(Date = dates[oos_idx], Lambda = lambda_oos)
p_lambda <- ggplot(lambda_df, aes(x = Date, y = log(Lambda))) +
  geom_line(color = "blue", linewidth = 0.8) +
  theme_classic() +
  labs(title = "Evolution of chosen λ (log scale) – ENet",
       x = "Time", y = "log(λ)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_lambda)

# α evolution
alpha_df <- data.frame(Date = dates[oos_idx], Alpha = alpha_oos)
p_alpha <- ggplot(alpha_df, aes(x = Date, y = Alpha)) +
  geom_line(color = "darkgreen", linewidth = 0.8) +
  theme_classic() +
  labs(title = "Evolution of chosen α – ENet",
       x = "Time", y = "α") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_alpha)

# ----------------------------------------------------------------------
# 14) Portfolio exercise using futures returns
# ----------------------------------------------------------------------
gamma <- 3
w_lower <- -1.5
w_upper <- 1.5
vol_window <- 60

asset_raw <- df$futures_ret
rf_raw <- df[[rf_col]]
rf_month <- (1 + rf_raw / 100)^(1/12) - 1
asset_excess <- asset_raw - rf_month

mu_enet <- enet_oos
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

valid <- is.finite(mu_enet) & is.finite(mu_bench) & is.finite(asset_excess_oos) &
  is.finite(rf_oos) & is.finite(sigma2)
mu_enet_p <- mu_enet[valid]
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

port_enet <- portfolio_metrics(mu_enet_p, asset_excess_p, rf_oos_p, sigma2_p)
port_bench <- portfolio_metrics(mu_bench_p, asset_excess_p, rf_oos_p, sigma2_p)
cer_gain_annual_pct <- 1200 * (port_enet$cer_monthly - port_bench$cer_monthly)

cat("\n--- Portfolio exercise (futures returns) ---\n")
cat("Annualized CER gain (ENet - Historical mean) [%] :", round(cer_gain_annual_pct, 4), "\n")
cat("Annualized Sharpe ratio (ENet)                   :", round(port_enet$sharpe_annual, 4), "\n")
cat("Annualized Sharpe ratio (Historical mean)         :", round(port_bench$sharpe_annual, 4), "\n")

wealth_df <- data.frame(Date = dates_p,
                        ENet = port_enet$wealth,
                        HistoricalMean = port_bench$wealth)
wealth_long <- melt(wealth_df, id.vars = "Date", variable.name = "Strategy", value.name = "Wealth")
p_wealth <- ggplot(wealth_long, aes(Date, Wealth, color = Strategy)) +
  geom_line(linewidth = 0.8) +
  theme_classic() +
  labs(title = "Cumulative wealth: Elastic Net vs Historical Mean",
       x = "Time", y = "Wealth (initial = 1)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_wealth)


