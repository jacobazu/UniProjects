

library(dplyr)
library(tidyr)
library(ggplot2)
library(forecast) 
library(rugarch)

# ----------------------------------------------------------------------
# 1. USER SETTINGS
# ----------------------------------------------------------------------
DATA_FILE  <- "ZhangData32_3.0.csv"
DATE_VAR   <- "date"
TARGET_VAR <- "target"        
RF_COL     <- "TB3MS"        
FUTURES_COL <- "Price_Futures"

# Split dates
TRAIN_END <- as.Date("1999-12-01")
VAL_END   <- as.Date("2007-12-01")
TEST_END  <- as.Date("2025-08-01")

# Rolling window length for predictor selection
WINDOW_SIZE <- 120

# Grid for tuning
pval_grid   <- c(0.01, 0.05, 0.10)
delta1_grid <- c(1, 1.5, 2)
delta2_grid <- c(1.5, 2, 2.5)
cv_scale <- 1

# Portfolio parameters
GAMMA <- 3
W_LOWER <- -1.5
W_UPPER <- 1.5
VOL_WINDOW <- 60   # months for rolling volatility

# ----------------------------------------------------------------------
# 2. PREDICTORS 
# ----------------------------------------------------------------------
tech_vars <- c("MA_1_9","MA_1_12","MA_2_9","MA_2_12","MA_3_9","MA_3_12",
               "MOM_1","MOM_2","MOM_3","MOM_6","MOM_9","MOM_12",
               "VOL_1_9","VOL_1_12","VOL_2_9","VOL_2_12","VOL_3_9","VOL_3_12")
macro_vars <- c("TB3MS","GS10","infl_m","SVOL","epu","kilian",
                "prod_growth","stocks_growth","imports_growth","m2_growth",
                "ip_growth","unemp_diff","cfnai","TCU","BAA","AAA")
predictors <- c(tech_vars, macro_vars)

# ----------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ----------------------------------------------------------------------
cp_threshold <- function(pval, delta, n, cscale = 1) {
  if (n <= 0) stop("n must be positive.")
  qnorm(1 - pval / (2 * cscale * (n ^ delta)))
}

safe_lm_tstat <- function(y, x, Xsel = NULL, Xcond = NULL) {
  y <- as.numeric(y); x <- as.numeric(x); n <- length(y)
  if (length(x) != n) stop("length mismatch")
  if (is.null(Xsel)) Xsel <- matrix(nrow = n, ncol = 0)
  else Xsel <- as.matrix(Xsel)
  if (is.null(Xcond)) Xcond <- matrix(nrow = n, ncol = 0)
  else Xcond <- as.matrix(Xcond)
  
  df <- data.frame(y = y, x = x)
  if (ncol(Xsel) > 0) {
    ctrl1 <- as.data.frame(Xsel); names(ctrl1) <- paste0("s", seq_len(ncol(ctrl1)))
    df <- cbind(df, ctrl1)
  }
  if (ncol(Xcond) > 0) {
    ctrl2 <- as.data.frame(Xcond); names(ctrl2) <- paste0("z", seq_len(ncol(ctrl2)))
    df <- cbind(df, ctrl2)
  }
  fit <- tryCatch(lm(y ~ ., data = df), error = function(e) NULL)
  if (is.null(fit)) return(-Inf)
  sm <- summary(fit)
  if (!("x" %in% rownames(sm$coefficients))) return(-Inf)
  abs(sm$coefficients["x", "t value"])
}

ocmt_select <- function(y, X, Xcond = NULL, pval = 0.05, delta1 = 1, delta2 = 2,
                        cscale = 1, max_stages = 50) {
  if (!is.matrix(X)) X <- as.matrix(X)
  n <- nrow(X); p <- ncol(X)
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

# Clark-West test for nested models
cw_test <- function(actual, bench, model) {
  e0 <- actual - bench
  e1 <- actual - model
  d <- e0^2 - (e1^2 - (bench - model)^2)
  n <- length(d)
  t_stat <- mean(d) / (sd(d) / sqrt(n))
  p_val <- 1 - pnorm(t_stat)
  c(t_stat = unname(t_stat), p_value = unname(p_val))
}

# Success ratio (directional accuracy on changes)
success_ratio <- function(actual, forecast) {
  actual_diff <- diff(actual)
  forecast_diff <- diff(forecast)
  if (length(forecast_diff) == length(actual_diff)) {
    correct <- sign(actual_diff) == sign(forecast_diff)
    return(mean(correct, na.rm = TRUE))
  } else {
    return(mean(sign(actual) == sign(forecast), na.rm = TRUE))
  }
}

# ----------------------------------------------------------------------
# 4. RECURSIVE FORECAST FUNCTION
# ----------------------------------------------------------------------
run_ocmt <- function(y_all, X_all, dates, start_idx, end_idx,
                     pval, delta1, delta2, window_size = 120, Xcond_all = NULL) {
  n_eval <- end_idx - start_idx + 1
  p <- ncol(X_all)
  f_model <- numeric(n_eval)
  f_bench <- numeric(n_eval)
  y_oos   <- numeric(n_eval)
  sel_mat <- matrix(0L, nrow = n_eval, ncol = p)
  colnames(sel_mat) <- colnames(X_all)
  
  for (ii in seq_len(n_eval)) {
    t <- start_idx + ii - 1
    

    sel_start <- max(1, t - window_size)
    sel_end <- t - 1
    sel_idx <- sel_start:sel_end
    y_sel <- y_all[sel_idx]
    X_sel <- X_all[sel_idx, , drop = FALSE]
    

    bench_start <- 1
    bench_end <- t - 1
    bench_idx <- bench_start:bench_end
    y_bench <- y_all[bench_idx]
    
    Xcond_sel <- if (!is.null(Xcond_all)) Xcond_all[sel_idx, , drop = FALSE] else NULL
    
    sel <- ocmt_select(y_sel, X_sel, Xcond_sel, pval, delta1, delta2, cv_scale)
    
    f0 <- mean(y_bench, na.rm = TRUE)            
    f1 <- forecast_ocmt(y_sel, X_sel, X_all[t, , drop = FALSE], sel$selected, Xcond_sel, NULL)
    
    y_oos[ii] <- y_all[t]
    f_bench[ii] <- f0
    f_model[ii] <- f1
    sel_mat[ii, ] <- as.integer(sel$selected)
  }
  list(dates = dates[start_idx:end_idx], y = y_oos, f_model = f_model,
       f_bench = f_bench, selection_matrix = sel_mat)
}

# ----------------------------------------------------------------------
# 5. LOAD AND PREPARE DATA
# ----------------------------------------------------------------------
data <- read.csv(DATA_FILE)
data[[DATE_VAR]] <- as.Date(data[[DATE_VAR]])

# Keep predictors, target, futures, risk‑free rate
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

# Compute futures simple returns 
df <- df %>%
  mutate(futures_ret = (lead(!!sym(FUTURES_COL)) - !!sym(FUTURES_COL)) / !!sym(FUTURES_COL)) %>%
  filter(!is.na(futures_ret))

# Update after removing last row
y_all <- y_all[1:nrow(df)]
dates <- dates[1:nrow(df)]
X_all <- X_all[1:nrow(df), , drop = FALSE]
rf_all <- df[[RF_COL]]
futures_ret_all <- df$futures_ret

# Split indices
train_end_idx <- max(which(dates <= TRAIN_END))
val_end_idx   <- max(which(dates <= VAL_END))
test_end_idx  <- max(which(dates <= TEST_END))

if (!is.finite(train_end_idx)) stop("TRAIN_END not found")
if (!is.finite(val_end_idx)) stop("VAL_END not found")
if (!is.finite(test_end_idx)) stop("TEST_END not found")

cat("Training   :", as.character(min(dates[1:train_end_idx])), "to", as.character(max(dates[1:train_end_idx])), "\n")
cat("Validation :", as.character(dates[train_end_idx+1]), "to", as.character(dates[val_end_idx]), "\n")
cat("Test       :", as.character(dates[val_end_idx+1]), "to", as.character(dates[test_end_idx]), "\n")

# ----------------------------------------------------------------------
# 6. GRID SEARCH ON VALIDATION 
# ----------------------------------------------------------------------
grid <- expand.grid(pval = pval_grid, delta1 = delta1_grid, delta2 = delta2_grid,
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid$rel_rmsfe <- NA_real_

for (g in seq_len(nrow(grid))) {
  res_val <- run_ocmt(y_all, X_all, dates,
                      start_idx = train_end_idx + 1, end_idx = val_end_idx,
                      pval = grid$pval[g], delta1 = grid$delta1[g], delta2 = grid$delta2[g],
                      window_size = WINDOW_SIZE)
  rmsfe_mod <- sqrt(mean((res_val$y - res_val$f_model)^2))
  rmsfe_ben <- sqrt(mean((res_val$y - res_val$f_bench)^2))
  grid$rel_rmsfe[g] <- rmsfe_mod / rmsfe_ben
  cat("grid", g, "of", nrow(grid), "done; rel RMSFE =", grid$rel_rmsfe[g], "\n")
}

best <- grid[which.min(grid$rel_rmsfe), ]
cat("\nBest hyperparameters:\n")
print(best)

# ----------------------------------------------------------------------
# 7. FINAL TEST EVALUATION
# ----------------------------------------------------------------------
final <- run_ocmt(y_all, X_all, dates,
                  start_idx = val_end_idx + 1, end_idx = test_end_idx,
                  pval = best$pval, delta1 = best$delta1, delta2 = best$delta2,
                  window_size = WINDOW_SIZE)

actual <- final$y           # log returns (percent)
pred_model <- final$f_model # log return forecasts (percent)
pred_bench <- final$f_bench # historical mean log return (percent)
test_dates <- final$dates

# Point forecast metrics
mspe_mod <- mean((actual - pred_model)^2)
mspe_ben <- mean((actual - pred_bench)^2)
oos_r2 <- 1 - mspe_mod / mspe_ben
rmsfe_mod <- sqrt(mspe_mod)
rmsfe_ben <- sqrt(mspe_ben)

# Directional metrics
sr <- success_ratio(actual, pred_model)
cw <- cw_test(actual, pred_bench, pred_model)
pt <- DACTest(diff(actual), diff(pred_model), test = "PT")

cat("\n========== FINAL TEST RESULTS ==========\n")
cat("RMSFE (OCMT):", rmsfe_mod, "\n")
cat("RMSFE (Benchmark):", rmsfe_ben, "\n")
cat("Relative RMSFE:", rmsfe_mod / rmsfe_ben, "\n")
cat("OOS R²:", oos_r2, "\n")
cat("Success ratio (on changes):", sr, "\n")
cat("Clark-West test: stat =", cw["t_stat"], "  p-value =", cw["p_value"], "\n")
cat("Pesaran-Timmermann test (on changes):\n")
print(pt)

# ----------------------------------------------------------------------
# 8. PORTFOLIO EXERCISE 
# ----------------------------------------------------------------------
# pred_model and pred_bench are already simple return forecasts (fractional)

# Risk‑free rate: annual percent to monthly simple
rf_month <- (1 + rf_all / 100)^(1/12) - 1

# Excess returns of futures (simple)
asset_excess <- futures_ret_all - rf_month

# Test period indices
test_rows <- (val_end_idx + 1):test_end_idx

# Rolling volatility (using past excess returns only)
sigma2 <- rep(NA_real_, length(test_rows))
for (t in seq_along(test_rows)) {
  global_t <- test_rows[t]
  if (global_t <= VOL_WINDOW) next
  hist_idx <- (global_t - VOL_WINDOW):(global_t - 1)
  sigma2[t] <- var(asset_excess[hist_idx], na.rm = TRUE)
}
# Forward fill first missing
first_sigma <- which(!is.na(sigma2))[1]
if (!is.na(first_sigma) && first_sigma > 1) sigma2[1:(first_sigma-1)] <- sigma2[first_sigma]

# Align all series – use pred_model and pred_bench directly
valid <- is.finite(pred_model) & is.finite(pred_bench) &
  is.finite(asset_excess[test_rows]) & is.finite(rf_month[test_rows]) & is.finite(sigma2)
pred_mod_p <- pred_model[valid]
pred_ben_p <- pred_bench[valid]
asset_excess_p <- asset_excess[test_rows][valid]
rf_p <- rf_month[test_rows][valid]
sigma2_p <- sigma2[valid]

# Portfolio weights
w_model <- (1 / GAMMA) * (pred_mod_p / sigma2_p)
w_bench <- (1 / GAMMA) * (pred_ben_p / sigma2_p)
w_model <- pmax(pmin(w_model, W_UPPER), W_LOWER)
w_bench <- pmax(pmin(w_bench, W_UPPER), W_LOWER)

# Portfolio returns
port_excess_model <- w_model * asset_excess_p
port_excess_bench <- w_bench * asset_excess_p
port_total_model <- rf_p + port_excess_model
port_total_bench <- rf_p + port_excess_bench

# CER  and Sharpe
cer_model <- mean(port_total_model) - 0.5 * GAMMA * var(port_total_model)
cer_bench <- mean(port_total_bench) - 0.5 * GAMMA * var(port_total_bench)
cer_gain_annual <- 1200 * (cer_model - cer_bench)

sharpe_model <- sqrt(12) * mean(port_excess_model) / sd(port_excess_model)
sharpe_bench <- sqrt(12) * mean(port_excess_bench) / sd(port_excess_bench)

cat("\n========== PORTFOLIO PERFORMANCE ==========\n")
cat("Annualized CER gain (OCMT - Historical mean) [%]:", round(cer_gain_annual, 4), "\n")
cat("Annualized Sharpe ratio (OCMT):", round(sharpe_model, 4), "\n")
cat("Annualized Sharpe ratio (Historical mean):", round(sharpe_bench, 4), "\n")

# Cumulative wealth plot
wealth_model <- cumprod(1 + port_total_model)
wealth_bench <- cumprod(1 + port_total_bench)
wealth_df <- data.frame(Date = dates[test_rows][valid],
                        OCMT = wealth_model,
                        HistoricalMean = wealth_bench)
wealth_long <- pivot_longer(wealth_df, cols = c("OCMT", "HistoricalMean"),
                            names_to = "Strategy", values_to = "Wealth")
p_wealth <- ggplot(wealth_long, aes(x = Date, y = Wealth, color = Strategy)) +
  geom_line(linewidth = 0.8) +
  theme_classic() +
  labs(title = "Cumulative Wealth: OCMT vs Historical Mean (simple returns)",
       x = "Time", y = "Wealth (initial = 1)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_wealth)
# ----------------------------------------------------------------------
# 9. CSPE DIFFERENCE PLOT
# ----------------------------------------------------------------------
cspe_diff <- cumsum((actual - pred_bench)^2 - (actual - pred_model)^2)
cspe_df <- data.frame(Date = test_dates, CSPE_Diff = cspe_diff)
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
# 10. MODEL SIZE OVER TIME
# ----------------------------------------------------------------------
model_size <- rowSums(final$selection_matrix)
size_df <- data.frame(Date = test_dates, ModelSize = model_size)
p_size <- ggplot(size_df, aes(x = Date, y = ModelSize)) +
  geom_line(color = "darkred", linewidth = 0.8) +
  geom_point(size = 0.8, color = "darkred") +
  theme_classic() +
  labs(title = "OCMT model size over time (rolling selection)",
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
# 11. VARIABLE SELECTION FREQUENCIES AND HEATMAP
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
    labs(title = "OCMT Selection Heatmap (rolling selection)", x = "Time", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 7), panel.grid = element_blank())
  print(p_heat)
}
