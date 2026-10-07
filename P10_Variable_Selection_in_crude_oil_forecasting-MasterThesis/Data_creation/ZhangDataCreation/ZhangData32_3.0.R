# constructing zhang data 

library(dplyr)
library(lubridate) 
library(tidyr)
library(quantmod)
library(eia)
library(readxl)
library(readr)
library(stringr)
library(TTR)
library(openxlsx)
library(tidyquant)   # for tq_get



# Helper: read FRED CSV and aggregate to monthly

read_fred_csv <- function(file_path, series_id) {
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  # First column is date (observation_date), second column is series value
  date_col <- names(df)[1]
  value_col <- names(df)[2]
  
  df <- df %>%
    mutate(
      date = as.Date(get(date_col)),
      value = as.numeric(get(value_col))
    ) %>%
    filter(!is.na(value)) %>%
    mutate(date = floor_date(date, "month")) %>%
    group_by(date) %>%
    summarise(value = last(value), .groups = "drop") %>%
    rename(!!series_id := value)
  
  return(df)
}


# 0. constructing futures data set

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


# 1. oil data (EIA spot)

oil_raw <- read_excel("RWTCm.xls")
oil_raw <- oil_raw %>%
  mutate(date = as.Date(Date, format = "%Y-%m-%d"))

oil_monthly <- oil_raw %>%
  rename(oil_price = oil_price_wti) %>%    # adjust column name as needed
  select(date, oil_price) %>%
  arrange(date)


# 2. align futures and spot to monthly (first day)

futures_monthly <- raw %>%
  mutate(
    yr = year(Date),
    dy = day(Date),      # currently stored as "month"
    mo = month(Date),    # currently stored as "day"
    date = as.Date(paste0(yr, "-", dy, "-", mo))
  ) %>%
  select(-yr, -mo, -dy)

futures_monthly <- futures_monthly %>%
  mutate(ym = floor_date(date, "month"))

spot_monthly <- oil_monthly %>%
  mutate(ym = floor_date(date, "month"))

prices <- full_join(
  futures_monthly %>% select(ym, Price_Futures, Vol),
  spot_monthly %>% select(ym, oil_price),
  by = "ym"
) %>%
  arrange(ym)


# 3. S&P 500 variance (from daily prices) 

sp500_daily_data <- tq_get("^GSPC",
                           from = "1980-01-01",
                           to = "2026-12-31",
                           get = "stock.prices")

SP500_monthly_SVOL <- sp500_daily_data %>%
  arrange(date) %>%
  mutate(
    daily_ret = (adjusted / lag(adjusted)) - 1,
    sq_ret = daily_ret^2
  ) %>%
  group_by(date = floor_date(date, "month")) %>%
  summarise(SVOL = sum(sq_ret, na.rm = TRUE)) %>%
  mutate(SVOL = lag(SVOL))   # lag so that SVOL in month t uses data up to t-1


# 4. macro variables from local CSV files

setwd("G:/Mit drev/AAU/10 SEM (speciale)/Actual/ZhangDataCreation/fredr datasets")
macro_files <- list(
  TB3MS      = "TB3MS.csv",
  GS10       = "GS10.csv",
  BAA        = "BAA.csv",
  AAA        = "AAA.csv",
  CPIAUCSL   = "CPIAUCSL.csv",
  INDPRO     = "INDPRO.csv",
  UNRATE     = "UNRATE.csv",
  M2SL       = "M2SL.csv",
  TCU        = "TCU.csv",
  IGREA      = "IGREA.csv",
  USEPUINDXM = "USEPUINDXM.csv",
  CFNAI      = "CFNAI.csv"
)

macro_list <- lapply(names(macro_files), function(id) {
  read_fred_csv(macro_files[[id]], id)
})

macro_monthly <- Reduce(function(x, y) full_join(x, y, by = "date"), macro_list) %>%
  arrange(date)

# Merge S&P 500 variance
macro_monthly <- macro_monthly %>%
  left_join(SP500_monthly_SVOL, by = "date")


# 5. EIA variables (oil production, imports, stocks)

setwd("G:/Mit drev/AAU/10 SEM (speciale)/Actual/ZhangDataCreation")
prod_path   <- "MCRFPUS2m.xls"
imports_path <- "MCRIMUS2m.xls"
stocks_path  <- "monthly_endingstocks.xlsx"

prod <- read_excel(prod_path) %>%
  rename(
    date = date,
    production = `production of oil (thousund of barrels per day)`
  ) %>%
  mutate(date = as.Date(date)) %>%
  group_by(date = floor_date(date, "month")) %>%
  summarise(production = last(production), .groups = "drop")

imports <- read_excel(imports_path) %>%
  rename(
    date = date,
    imports = `import of crude oil (thousand of barrels per day)`
  ) %>%
  mutate(date = as.Date(date)) %>%
  group_by(date = floor_date(date, "month")) %>%
  summarise(imports = last(imports), .groups = "drop")

stocks <- read_excel(stocks_path) %>%
  rename(
    date = date,
    stocks = endingstocks
  ) %>%
  mutate(date = as.Date(date)) %>%
  group_by(date = floor_date(date, "month")) %>%
  summarise(stocks = last(stocks), .groups = "drop")

macro_monthly <- macro_monthly %>%
  left_join(prod, by = "date") %>%
  left_join(imports, by = "date") %>%
  left_join(stocks, by = "date")


# 6. transform macro variables (growth rates, lags)

macro_monthly <- macro_monthly %>%
  arrange(date) %>%
  mutate(
    infl_m = 100 * (log(CPIAUCSL) - log(lag(CPIAUCSL))),
    m2_growth = 100 * (log(M2SL) - log(lag(M2SL))),
    ip_growth = 100 * (log(INDPRO) - log(lag(INDPRO))),
    unemp_diff = UNRATE - lag(UNRATE),
    prod_growth = 100 * (log(production) - log(lag(production))),
    stocks_growth = 100 * (log(stocks) - log(lag(stocks))),
    imports_growth = 100 * (log(imports) - log(lag(imports))),
    kilian = IGREA,
    epu = USEPUINDXM,
    cfnai = CFNAI
  )

# Lag all except CPI (CPIAUCSL stays as is for deflation)
macro_monthly <- macro_monthly %>%
  arrange(date) %>%
  mutate(across(-c(date, CPIAUCSL), ~ lag(.x, 1)))


# 7. merge oil data with macro data

prices <- prices %>% rename(date = ym)

full_data <- prices %>%
  left_join(macro_monthly, by = "date") %>%
  arrange(date) %>%
  mutate(
    real_price = oil_price / CPIAUCSL * 100,
    spot_ret_simple = (real_price - lag(real_price)) / lag(real_price),
    target = lead(spot_ret_simple)
  )


# 8. technical indicators (monthly) from futures data

MA_signal <- function(price, short, long) {
  sma <- SMA(price, n = short)
  lma <- SMA(price, n = long)
  ifelse(sma >= lma, 1, 0)
}

MOM_signal <- function(price, m) {
  ifelse(price >= lag(price, m), 1, 0)
}

VOL_signal <- function(price, vol, short, long) {
  obv <- OBV(price, vol)
  sma_obv <- SMA(obv, n = short)
  lma_obv <- SMA(obv, n = long)
  ifelse(sma_obv >= lma_obv, 1, 0)
}

tech_monthly <- full_data %>%
  select(date, Price_Futures, Vol) %>%
  arrange(date) %>%
  mutate(
    MA_1_9  = MA_signal(Price_Futures, 1, 9),
    MA_1_12 = MA_signal(Price_Futures, 1, 12),
    MA_2_9  = MA_signal(Price_Futures, 2, 9),
    MA_2_12 = MA_signal(Price_Futures, 2, 12),
    MA_3_9  = MA_signal(Price_Futures, 3, 9),
    MA_3_12 = MA_signal(Price_Futures, 3, 12),
    
    MOM_1  = MOM_signal(Price_Futures, 1),
    MOM_2  = MOM_signal(Price_Futures, 2),
    MOM_3  = MOM_signal(Price_Futures, 3),
    MOM_6  = MOM_signal(Price_Futures, 6),
    MOM_9  = MOM_signal(Price_Futures, 9),
    MOM_12 = MOM_signal(Price_Futures, 12),
    
    VOL_1_9  = VOL_signal(Price_Futures, Vol, 1, 9),
    VOL_1_12 = VOL_signal(Price_Futures, Vol, 1, 12),
    VOL_2_9  = VOL_signal(Price_Futures, Vol, 2, 9),
    VOL_2_12 = VOL_signal(Price_Futures, Vol, 2, 12),
    VOL_3_9  = VOL_signal(Price_Futures, Vol, 3, 9),
    VOL_3_12 = VOL_signal(Price_Futures, Vol, 3, 12)
  ) %>%
  select(date, starts_with("MA_"), starts_with("MOM_"), starts_with("VOL_"))

full_data <- full_data %>%
  left_join(tech_monthly, by = "date") %>%
  arrange(date)

head(full_data)
write.csv(full_data, "ZhangData32_3.0.csv", row.names = FALSE)
