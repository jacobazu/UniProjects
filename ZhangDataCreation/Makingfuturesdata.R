
# constructing futures data


# 1. reading futures csv (the sourceis investing.com)
raw <- read.csv("futures1983-2026.csv", stringsAsFactors = FALSE)

# 2. convert to dd/mm/yyyy
raw$Date <- as.Date(raw$Date, format = "%d/%m/%Y")

# 3. convert price to numeric
raw$Price <- as.numeric(raw$Price)

# 4. convert Vol to actual numbers
convert_vol <- function(x) {
  x <- gsub(",", "", x)       
  x <- gsub("M", "*1e6", x)   
  x <- gsub("K", "*1e3", x)   
  sapply(x, function(y) eval(parse(text = y)))  
}

raw$Vol <- convert_vol(raw$Vol)

# 5. keep only relevant columns
raw <- raw[, c("Date", "Price", "Vol")]

# 6. rename Price to Price_Futures
names(raw)[names(raw) == "Price"] <- "Price_Futures"

# 7. checking
head(raw)
