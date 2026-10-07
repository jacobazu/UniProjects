
library(dplyr)
library(tidyr)
library(ggplot2)
library(rugarch)

# ---------------------------
# 1. USER SETTINGS
# ---------------------------
DATA_FILE  <- "ZhangData32_3.0.csv"
DATE_VAR   <- "date"
TARGET_VAR <- "target"
RF_COL     <- "TB3MS"         
FUTURES_COL <- "Price_Futures"  

TRAIN_END <- as.Date("1999-12-01")
TEST_END  <- as.Date("2025-08-01")

# Fixed hyperparameters
PVAL   <- 0.01
DELTA1 <- 1
DELTA2 <- 2

# Portfolio parameters
GAMMA <- 3
W_LOWER <- -1.5
W_UPPER <- 1.5
VOL_WINDOW <- 60   

# ---------------------------
# 2. PREDICTORS
# ---------------------------
tech_vars <- c("MA_1_9","MA_1_12","MA_2_9","MA_2_12","MA_3_9","MA_3_12",
               "MOM_1","MOM_2","MOM_3","MOM_6","MOM_9","MOM_12",
               "VOL_1_9","VOL_1_12","VOL_2_9","VOL_2_12","VOL_3_9","VOL_3_12")

macro_vars <- c("TB3MS","GS10","infl_m","SVOL","epu","kilian",
                "prod_growth","stocks_growth","imports_growth","m2_growth",
                "ip_growth","unemp_diff","cfnai","TCU","BAA","AAA")

predictors <- c(tech_vars, macro_vars)

# ---------------------------
# 3. HELPER FUNCTIONS
# ---------------------------
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

clark_west_test <- function(y, f_bench, f_model) {
  d <- (y - f_bench)^2 - ((y - f_model)^2 - (f_bench - f_model)^2)
  n <- length(d)
  t_stat <- mean(d) / (sd(d) / sqrt(n))
  p_val  <- 1 - pnorm(t_stat)
  list(stat = t_stat, p_value = p_val)
}

pesaran_timmermann_test <- function(y, f) {
  actual <- as.integer(y > 0)
  pred   <- as.integer(f > 0)
  
  p_hat <- mean(actual)
  q_hat <- mean(pred)
  hit   <- mean(actual == pred)
  
  hit_0 <- p_hat * q_hat + (1 - p_hat) * (1 - q_hat)
  var_0  <- (p_hat * (1 - p_hat) * q_hat * (1 - q_hat)) / length(y)
  
  if (!is.finite(var_0) || var_0 <= 0) {
    return(list(stat = NA_real_, p_value = NA_real_, hit_rate = hit))
  }
  
  stat <- (hit - hit_0) / sqrt(var_0)
  p_val <- 1 - pnorm(stat)
  
  list(stat = stat, p_value = p_val, hit_rate = hit)
}

abs_t_first <- function(y, x, Xsel = NULL) {
  x <- as.numeric(x)
  n <- length(y)
  
  if (is.null(Xsel)) {
    Xsel <- matrix(nrow = n, ncol = 0)
  } else {
    Xsel <- as.matrix(Xsel)
    if (ncol(Xsel) == 0) {
      Xsel <- matrix(nrow = n, ncol = 0)
    }
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

# ---------------------------
# 4. BMT SELECTION
# ---------------------------
boosting_glm <- function(y, X, pval, delta1, delta2) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (nrow(X) != length(y)) stop("length(y) must equal nrow(X).")
  
  N <- ncol(X)
  if (N < 1) stop("X must have at least one column.")
  
  p1 <- pval / (N^(delta1 - 1))
  p2 <- pval / (N^(delta2 - 1))
  
  t1 <- qnorm(1 - p1 / (2 * N))
  t2 <- qnorm(1 - p2 / (2 * N))
  
  ind <- rep(FALSE, N)
  
  # initial selection
  ts <- rep(-Inf, N)
  for (i in seq_len(N)) {
    ts[i] <- abs_t_first(y = y, x = X[, i, drop = TRUE], Xsel = NULL)
  }
  
  i_max <- which.max(ts)
  if (is.finite(ts[i_max]) && ts[i_max] > t1) {
    ind[i_max] <- TRUE
  }
  
  # boosting loop
  repeat {
    ts <- rep(-Inf, N)
    
    for (j in seq_len(N)) {
      if (!ind[j]) {
        ts[j] <- abs_t_first(
          y = y,
          x = X[, j, drop = TRUE],
          Xsel = X[, ind, drop = FALSE]
        )
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

# ---------------------------
# 5. ONE-STEP-AHEAD FORECAST
# ---------------------------
forecast_one_step <- function(y_train, X_train, X_test, selected) {
  if (!any(selected)) {
    return(mean(y_train, na.rm = TRUE))
  }
  
  Xsel_train <- as.data.frame(X_train[, selected, drop = FALSE])
  Xsel_test  <- as.data.frame(t(X_test[selected]))
  
  if (ncol(Xsel_train) == 0) {
    return(mean(y_train, na.rm = TRUE))
  }
  
  names(Xsel_train) <- names(Xsel_test) <- colnames(X_train)[selected]
  
  fit <- tryCatch({
    lm(y_train ~ ., data = data.frame(y_train = y_train, Xsel_train))
  }, error = function(e) NULL)
  
  if (is.null(fit)) {
    return(mean(y_train, na.rm = TRUE))
  }
  
  pred <- tryCatch({
    as.numeric(predict(fit, newdata = Xsel_test))
  }, error = function(e) NA_real_)
  
  if (!is.finite(pred)) mean(y_train, na.rm = TRUE) else pred
}

# ---------------------------
# 6. RECURSIVE EVALUATION (expanding window, no validation)
# ---------------------------
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
    
    selected <- boosting_glm(
      y = y_train,
      X = X_train_sc,
      pval = pval,
      delta1 = d1,
      delta2 = d2
    )
    
    f0 <- mean(y_train, na.rm = TRUE)
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

# ---------------------------
# 7. LOAD AND PREPARE DATA
# ---------------------------
data <- read.csv(DATA_FILE)
data[[DATE_VAR]] <- as.Date(data[[DATE_VAR]])

keep_cols <- c(DATE_VAR, TARGET_VAR, RF_COL, FUTURES_COL, predictors)
df <- data %>%
  select(all_of(keep_cols)) %>%
  arrange(.data[[DATE_VAR]]) %>%
  filter(.data[[DATE_VAR]] >= as.Date("1986-02-01"),
         .data[[DATE_VAR]] <= TEST_END) %>%
  na.omit()

dates <- df[[DATE_VAR]]
y_all <- df[[TARGET_VAR]]
X_all <- as.matrix(df[, predictors, drop = FALSE])
colnames(X_all) <- predictors

# Futures returns 
df <- df %>%
  mutate(futures_ret = (lead(!!sym(FUTURES_COL)) - !!sym(FUTURES_COL)) / !!sym(FUTURES_COL)) %>%
  filter(!is.na(futures_ret))

# Update after removing last row
y_all <- y_all[1:nrow(df)]
dates <- dates[1:nrow(df)]
X_all <- X_all[1:nrow(df), , drop = FALSE]
rf_all <- df[[RF_COL]]
futures_ret_all <- df$futures_ret

# Define indices: training ends at TRAIN_END, test starts next month
train_end_idx <- max(which(dates <= TRAIN_END))
test_start_idx <- train_end_idx + 1
test_end_idx   <- max(which(dates <= TEST_END))

cat("Training period  :", as.character(min(dates[1:train_end_idx])), "to", as.character(max(dates[1:train_end_idx])), "\n")
cat("Test period      :", as.character(dates[test_start_idx]), "to", as.character(dates[test_end_idx]), "\n")
cat("Fixed hyperparameters: pval =", PVAL, ", delta1 =", DELTA1, ", delta2 =", DELTA2, "\n")

# ---------------------------
# 8. RUN BMT ON TEST PERIOD
# ---------------------------
final <- run_forecast(
  start_idx = test_start_idx,
  end_idx   = test_end_idx,
  y_all     = y_all,
  X_all     = X_all,
  pval      = PVAL,
  d1        = DELTA1,
  d2        = DELTA2
)

# ---------------------------
# 9. TEST METRICS
# ---------------------------
cw <- clark_west_test(final$y, final$f_bench, final$f_model)
pt <- pesaran_timmermann_test(final$y, final$f_model)

cat("\n========== FINAL TEST RESULTS ==========\n")
cat("OOS R²:", final$r2, "\n")
cat("Clark-West stat:", cw$stat, "  p-value:", cw$p_value, "\n")
cat("Pesaran-Timmermann stat:", pt$stat, "  p-value:", pt$p_value, "  hit rate:", pt$hit_rate, "\n")

# ---------------------------
# 10. MODEL SIZE OVER TIME
# ---------------------------
dates_test <- dates[test_start_idx:test_end_idx]
model_size <- rowSums(final$sel)
model_size_df <- data.frame(Date = dates_test, ModelSize = model_size)

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

# ---------------------------
# 11. CSPE DIFFERENCE PLOT
# ---------------------------
cspe_diff <- cumsum((final$y - final$f_bench)^2 - (final$y - final$f_model)^2)
cspe_df <- data.frame(Date = dates_test, CSPE_Diff = cspe_diff)

p_cspe <- ggplot(cspe_df, aes(x = Date, y = CSPE_Diff)) +
  geom_line(linewidth = 0.8, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred") +
  theme_classic() +
  labs(title = paste("Cumulative squared prediction error difference (BMT vs Historical mean)\n",
                     "pval =", PVAL, ", δ1 =", DELTA1, ", δ2 =", DELTA2),
       x = "Time", y = "CSPE (Benchmark - BMT)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_cspe)

# ---------------------------
# 12. SELECTION HEATMAP
# ---------------------------
heat_df <- as.data.frame(final$sel)
heat_df$date <- dates_test
long_df <- heat_df %>%
  pivot_longer(cols = all_of(predictors), names_to = "variable", values_to = "selected")

p_heat <- ggplot(long_df, aes(x = date, y = variable, fill = factor(selected))) +
  geom_tile() +
  scale_fill_manual(values = c("0" = "white", "1" = "blue"), name = "Selected") +
  labs(title = "BMT Selection Heatmap", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 7),
        panel.grid = element_blank())
print(p_heat)

# ---------------------------
# 13. VARIABLE IMPORTANCE
# ---------------------------
importance <- colMeans(final$sel)
importance_sorted <- sort(importance, decreasing = TRUE)
cat("\n--- Variable selection frequencies (%) ---\n")
print(round(100 * importance_sorted, 1))

# ---------------------------
# 14. PORTFOLIO EXERCISE
# ---------------------------
test_rows <- test_start_idx:test_end_idx

rf_month <- (1 + rf_all / 100)^(1/12) - 1
asset_excess <- futures_ret_all - rf_month

mu_bmt <- final$f_model
mu_bench <- final$f_bench

# Rolling volatility
sigma2 <- rep(NA_real_, length(asset_excess[test_rows]))
for (t in seq_along(test_rows)) {
  global_t <- test_rows[t]
  if (global_t <= VOL_WINDOW) next
  hist_idx <- (global_t - VOL_WINDOW):(global_t - 1)
  sigma2[t] <- var(asset_excess[hist_idx], na.rm = TRUE)
}
first_sigma <- which(!is.na(sigma2))[1]
if (!is.na(first_sigma) && first_sigma > 1) sigma2[1:(first_sigma-1)] <- sigma2[first_sigma]

# Align valid observations
valid <- is.finite(mu_bmt) & is.finite(mu_bench) & is.finite(asset_excess[test_rows]) &
  is.finite(rf_month[test_rows]) & is.finite(sigma2)
mu_bmt_p <- mu_bmt[valid]
mu_bench_p <- mu_bench[valid]
asset_excess_p <- asset_excess[test_rows][valid]
rf_p <- rf_month[test_rows][valid]
sigma2_p <- sigma2[valid]
dates_p <- dates[test_rows][valid]

# Portfolio weights
w_bmt <- (1 / GAMMA) * (mu_bmt_p / sigma2_p)
w_bench <- (1 / GAMMA) * (mu_bench_p / sigma2_p)
w_bmt <- pmax(pmin(w_bmt, W_UPPER), W_LOWER)
w_bench <- pmax(pmin(w_bench, W_UPPER), W_LOWER)

# Portfolio returns
port_excess_bmt <- w_bmt * asset_excess_p
port_excess_bench <- w_bench * asset_excess_p
port_total_bmt <- rf_p + port_excess_bmt
port_total_bench <- rf_p + port_excess_bench

# Cumulative wealth
wealth_bmt <- cumprod(1 + port_total_bmt)
wealth_bench <- cumprod(1 + port_total_bench)

# CER and Sharpe
cer_bmt <- mean(port_total_bmt) - 0.5 * GAMMA * var(port_total_bmt)
cer_bench <- mean(port_total_bench) - 0.5 * GAMMA * var(port_total_bench)
cer_gain_annual <- 1200 * (cer_bmt - cer_bench)

sharpe_bmt <- sqrt(12) * mean(port_excess_bmt) / sd(port_excess_bmt)
sharpe_bench <- sqrt(12) * mean(port_excess_bench) / sd(port_excess_bench)

cat("\n========== PORTFOLIO PERFORMANCE ==========\n")
cat("Annualized CER gain (BMT - Historical mean) [%]:", round(cer_gain_annual, 4), "\n")
cat("Annualized Sharpe ratio (BMT):", round(sharpe_bmt, 4), "\n")
cat("Annualized Sharpe ratio (Historical mean):", round(sharpe_bench, 4), "\n")

# Cumulative wealth plot
wealth_df <- data.frame(Date = dates_p,
                        BMT = wealth_bmt,
                        HistoricalMean = wealth_bench)
wealth_long <- pivot_longer(wealth_df, cols = c("BMT", "HistoricalMean"),
                            names_to = "Strategy", values_to = "Wealth")

p_wealth <- ggplot(wealth_long, aes(x = Date, y = Wealth, color = Strategy)) +
  geom_line(linewidth = 0.8) +
  theme_classic() +
  labs(title = "Cumulative Wealth: BMT vs Historical Mean",
       x = "Time", y = "Wealth (initial = 1)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_wealth)