# ============================================================================
# Complete dataset construction for crude oil volatility forecasting
#   - Target: LV = log(realized volatility) from daily WTI spot prices
#   - Predictors: 14 uncertainty variables + technical indicators 
#   - Additional uncertainty indices for exploratory analysis
# ============================================================================


library(tidyquant)
library(dplyr)
library(lubridate)
library(readxl)
library(tseries)
library(TTR)


# 1. Daily WTI spot prices -> LV 

wti_daily <- read.csv("datasets/DCOILWTICO.csv") %>%
  mutate(date = as.Date(date)) %>%
  arrange(date) %>%
  filter(!is.na(price), price > 0)

wti_daily <- wti_daily %>%
  mutate(ret = log(price / lag(price)),
         sq_ret = ret^2)

LV <- wti_daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(RV = sum(sq_ret, na.rm = TRUE), .groups = "drop") %>%
  mutate(LV = log(RV)) %>%
  rename(date = month)


# 2. Daily futures data for technical indicators

raw <- read.csv("futures1983-2026.csv", stringsAsFactors = FALSE)
raw$Date <- as.Date(raw$Date, format = "%d/%m/%Y")   
raw$Price <- as.numeric(raw$Price)

convert_vol <- function(x) {
  x <- gsub(",", "", x)
  x <- gsub("M", "*1e6", x)
  x <- gsub("K", "*1e3", x)
  sapply(x, function(y) eval(parse(text = y)))
}
raw$Vol <- convert_vol(raw$Vol)

raw <- raw[, c("Date", "Price", "Vol")]
names(raw)[names(raw) == "Price"] <- "Price_Futures"
raw <- raw[order(raw$Date), ]
rownames(raw) <- NULL

raw <- raw %>%
  mutate(
    yr = year(Date),
    dy = day(Date),      # currently stored as "month"
    mo = month(Date),    # currently stored as "day"
    date = as.Date(paste0(yr, "-", dy, "-", mo))
  ) %>%
  select(-yr, -mo, -dy)

raw <- raw %>%
  mutate(ym = floor_date(date, "month"))
raw <- raw %>% select(-Date)   
raw <- raw %>% select(-date)  
raw <- raw %>% rename(Date = ym)

raw <- raw %>% relocate(Date, .before = 1)


# 3. Construct technical indicators on monthly futures

raw <- raw %>%
  mutate(
    # Moving averages 
    MA_1_9  = as.integer(SMA(Price_Futures, 1) >= SMA(Price_Futures, 9)),
    MA_1_12 = as.integer(SMA(Price_Futures, 1) >= SMA(Price_Futures, 12)),
    MA_2_9  = as.integer(SMA(Price_Futures, 2) >= SMA(Price_Futures, 9)),
    MA_2_12 = as.integer(SMA(Price_Futures, 2) >= SMA(Price_Futures, 12)),
    MA_3_9  = as.integer(SMA(Price_Futures, 3) >= SMA(Price_Futures, 9)),
    MA_3_12 = as.integer(SMA(Price_Futures, 3) >= SMA(Price_Futures, 12)),
    
    # Momentum
    MOM_1  = as.integer(Price_Futures >= lag(Price_Futures, 1)),
    MOM_2  = as.integer(Price_Futures >= lag(Price_Futures, 2)),
    MOM_3  = as.integer(Price_Futures >= lag(Price_Futures, 3)),
    MOM_6  = as.integer(Price_Futures >= lag(Price_Futures, 6)),
    MOM_9  = as.integer(Price_Futures >= lag(Price_Futures, 9)),
    MOM_12 = as.integer(Price_Futures >= lag(Price_Futures, 12)),
    
    # OBV 
    OBV = OBV(Price_Futures, Vol),
    VOL_1_9  = as.integer(SMA(OBV, 1) >= SMA(OBV, 9)),
    VOL_1_12 = as.integer(SMA(OBV, 1) >= SMA(OBV, 12)),
    VOL_2_9  = as.integer(SMA(OBV, 2) >= SMA(OBV, 9)),
    VOL_2_12 = as.integer(SMA(OBV, 2) >= SMA(OBV, 12)),
    VOL_3_9  = as.integer(SMA(OBV, 3) >= SMA(OBV, 9)),
    VOL_3_12 = as.integer(SMA(OBV, 3) >= SMA(OBV, 12))
  )

# Remove initial rows with NAs 
raw_clean <- raw %>%
  filter(complete.cases(select(., starts_with("MA_"), starts_with("MOM_"))))

# Keep only the technical indicators 
tech_monthly <- raw_clean %>%
  select(Date, starts_with("MA_"), starts_with("MOM_"), starts_with("VOL_"))

# Rename Date to date for merging
tech_monthly <- tech_monthly %>% rename(date = Date)


# 4. Uncertainty variables (core 14 + additional)

# EPU
epu_raw <- read.csv("datasets/USEPUINDXM.csv")
epu <- epu_raw %>%
  select(date = month, EPU = USEPUINDXM) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

# GPR
gpr_raw <- read_excel("datasets/GPR.xlsx")
gpr <- gpr_raw %>%
  select(date = date, GPR = GPR) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

# EMV
emv_raw <- read_excel("datasets/EMV.xlsx")
emv <- emv_raw %>%
  select(date = date, EMV = EMV) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

# MPU
mpu_raw <- read_excel("datasets/MPU.xlsx")
mpu <- mpu_raw %>%
  select(date = date, MPU = MPU) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

# S&P 500 monthly variance and skewness
sp500 <- tq_get("^GSPC", get = "stock.prices", from = "1985-01-01")
sp500_monthly <- sp500 %>%
  mutate(ret = (adjusted / lag(adjusted)) - 1) %>%
  group_by(date = floor_date(date, "month")) %>%
  summarise(
    SVAR = sum(ret^2, na.rm = TRUE),
    RSK_SP = skewness(ret, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rename(date = date)

# RSK_Oil from daily WTI returns
rsk_oil <- wti_daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month) %>%
  summarise(RSK_Oil = skewness(ret, na.rm = TRUE), .groups = "drop") %>%
  rename(date = month)

# RA and EU (from Nancy Xu)
ra_raw <- read_excel("datasets/RA.xlsx")
ra_and_unc <- ra_raw %>%
  select(date = date, RA = RA, EU = UNC) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

# RU, MU, FU (Jurado-Ludvigson-Ng)
mu_ru_fu_raw <- read_excel("datasets/MURUFU.xlsx")
mu_ru_fu <- mu_ru_fu_raw %>%
  select(date = date, RU = RU, MU = MU, FU = FU) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

# VIX
vix <- tq_get("^VIX", get = "stock.prices", from = "1990-01-01") %>%
  group_by(date = floor_date(date, "month")) %>%
  summarise(VIX = last(adjusted))

# OVX
ovx <- tq_get("^OVX", get = "stock.prices", from = "1990-01-01") %>%
  group_by(date = floor_date(date, "month")) %>%
  summarise(OVX = last(adjusted))

# Additional (optional) uncertainty indices
cpu_raw <- read_excel("datasets/CPU.xlsx")
cpu <- cpu_raw %>%
  select(date = date, CPU = CPU) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

gepu_raw <- read_excel("datasets/GEPU.xlsx")
gepu <- gepu_raw %>%
  select(date = date, GEPU = GEPU) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

wui_raw <- read_excel("datasets/WUI.xlsx")
wui <- wui_raw %>%
  select(date = date, WUI = WUI) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)

idemv_raw <- read.csv("datasets/IDEMV.csv")
idemv_raw <- idemv_raw %>%
  mutate(date = make_date(year, month, day)) %>%
  arrange(date)
idemv <- idemv_raw %>%
  group_by(date = floor_date(date, "month")) %>%
  summarise(IDEMV = last(daily_infect_emv_index), .groups = "drop")

tpu_raw <- read_excel("datasets/TPU.xlsx")
tpu <- tpu_raw %>%
  select(date = date, TPU = TPU) %>%
  mutate(date = as.Date(date)) %>%
  arrange(date)


# 5. Helper: fix all dates to month start 

fix_dates <- function(df, col = "date") {
  if (!col %in% names(df)) return(df)
  x <- df[[col]]
  if (is.numeric(x)) x <- as.Date(x, origin = "1899-12-30")
  else if (inherits(x, "POSIXct")) x <- as.Date(x)
  else if (is.character(x)) x <- as.Date(x)
  x <- floor_date(x, "month")
  df[[col]] <- x
  df
}

LV <- fix_dates(LV)
epu <- fix_dates(epu)
gpr <- fix_dates(gpr)
emv <- fix_dates(emv)
mpu <- fix_dates(mpu)
sp500_monthly <- fix_dates(sp500_monthly)
rsk_oil <- fix_dates(rsk_oil)
ra_and_unc <- fix_dates(ra_and_unc)
mu_ru_fu <- fix_dates(mu_ru_fu)
vix <- fix_dates(vix)
ovx <- fix_dates(ovx)
cpu <- fix_dates(cpu)
gepu <- fix_dates(gepu)
wui <- fix_dates(wui)
idemv <- fix_dates(idemv)
tpu <- fix_dates(tpu)
tech_monthly <- fix_dates(tech_monthly)


# 6. Merge everything 

full_data <- LV %>%
  full_join(epu, by = "date") %>%
  full_join(gpr, by = "date") %>%
  full_join(emv, by = "date") %>%
  full_join(mpu, by = "date") %>%
  full_join(sp500_monthly, by = "date") %>%
  full_join(rsk_oil, by = "date") %>%
  full_join(ra_and_unc, by = "date") %>%
  full_join(mu_ru_fu, by = "date") %>%
  full_join(vix, by = "date") %>%
  full_join(ovx, by = "date") %>%
  full_join(cpu, by = "date") %>%
  full_join(gepu, by = "date") %>%
  full_join(wui, by = "date") %>%
  full_join(idemv, by = "date") %>%
  full_join(tpu, by = "date") %>%
  full_join(tech_monthly, by = "date") %>%
  arrange(date)


# 7. Apply transformations on the full dataset

full_data <- full_data %>%
  mutate(
    # First differences for non‑stationary series
    d_EPU   = EPU - lag(EPU),
    d_MPU   = MPU - lag(MPU),
    d_GEPU  = GEPU - lag(GEPU),
    d_TPU   = TPU - lag(TPU),
    d_CPU   = CPU - lag(CPU),
    d_RU    = RU - lag(RU),
    d_MU    = MU - lag(MU),
    d_FU    = FU - lag(FU),
    d_WUI   = WUI - lag(WUI),
    d_IDEMV = IDEMV - lag(IDEMV),
    
    # Levels (stationary or trend‑stationary)
    EMV     = EMV,
    SVAR    = SVAR,
    RSK_SP  = RSK_SP,
    RSK_Oil = RSK_Oil,
    GPR     = GPR,
    RA      = RA,
    EU      = EU,
    
    # Log levels
    log_VIX = log(VIX),
    log_OVX = log(OVX),
    
    # Technical indicators 
    MA_1_9  = MA_1_9,
    MA_1_12 = MA_1_12,
    MA_2_9  = MA_2_9,
    MA_2_12 = MA_2_12,
    MA_3_9  = MA_3_9,
    MA_3_12 = MA_3_12,
    MOM_1   = MOM_1,
    MOM_2   = MOM_2,
    MOM_3   = MOM_3,
    MOM_6   = MOM_6,
    MOM_9   = MOM_9,
    MOM_12  = MOM_12,
    VOL_1_9  = VOL_1_9,
    VOL_1_12 = VOL_1_12,
    VOL_2_9  = VOL_2_9,
    VOL_2_12 = VOL_2_12,
    VOL_3_9  = VOL_3_9,
    VOL_3_12 = VOL_3_12
  )


# 8. Lag all predictors by one month 

predictor_cols <- c("d_EPU", "d_MPU", "d_GEPU", "d_TPU", "d_CPU",
                    "d_RU", "d_MU", "d_FU", "d_WUI", "d_IDEMV",
                    "EMV", "SVAR", "RSK_SP", "RSK_Oil", "GPR", "RA", "EU",
                    "log_VIX", "log_OVX",
                    "MA_1_9", "MA_1_12", "MA_2_9", "MA_2_12", "MA_3_9", "MA_3_12",
                    "MOM_1", "MOM_2", "MOM_3", "MOM_6", "MOM_9", "MOM_12",
                    "VOL_1_9", "VOL_1_12", "VOL_2_9", "VOL_2_12", "VOL_3_9", "VOL_3_12")

full_data <- full_data %>%
  mutate(across(all_of(predictor_cols), ~ lag(.x, 1), .names = "{.col}_lag"))


# 9. Filter to 2008-2024

full_data_since_2008 <- full_data %>%
  filter(date >= as.Date("2008-02-01") & date <= as.Date("2024-02-01"))


# 10. Final modelling dataset 

modelling_data <- full_data_since_2008 %>%
  select(date, LV, ends_with("_lag")) %>%
  filter(complete.cases(.))

# Remove "_lag" suffix for simpler names
names(modelling_data) <- gsub("_lag$", "", names(modelling_data))


# 11. Sample summary

cat("Sample start:", as.character(min(modelling_data$date)), "\n")
cat("Sample end  :", as.character(max(modelling_data$date)), "\n")
cat("Observations:", nrow(modelling_data), "\n")


# 12. Stationarity check after transformations

transformed_vars <- predictor_cols
test_data <- full_data_since_2008 %>% filter(complete.cases(across(all_of(transformed_vars))))
stationarity_results <- data.frame(Variable = character(),
                                   ADF_p = numeric(),
                                   KPSS_p = numeric(),
                                   Stationary = logical())

for (var in transformed_vars) {
  series <- test_data[[var]]
  adf_p <- tryCatch(adf.test(series, alternative = "stationary")$p.value, error = function(e) NA)
  kpss_p <- tryCatch(kpss.test(series, null = "Level")$p.value, error = function(e) NA)
  stationary <- ifelse(!is.na(adf_p) & !is.na(kpss_p), adf_p < 0.05 & kpss_p > 0.05, FALSE)
  stationarity_results <- rbind(stationarity_results,
                                data.frame(Variable = var, ADF_p = adf_p,
                                           KPSS_p = kpss_p, Stationary = stationary))
}
cat("\n--- Stationarity check after transformations ---\n")
print(stationarity_results)


# 13. Save final dataset

write.csv(modelling_data, "VolatilityData.csv", row.names = FALSE)

