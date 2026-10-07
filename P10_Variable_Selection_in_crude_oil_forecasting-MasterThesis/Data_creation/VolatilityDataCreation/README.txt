# Data generation

R scripts that build the datasets for my thesis on forecasting crude oil
returns and volatility.

## What's here

- `futures_data.R` builds the monthly WTI futures series (price and volume).
  Both datasets use it for the technical indicators.
- `ZhangDataCreation/` has the script for the returns dataset.
- `VolatilityDataCreation/` has the script for the volatility dataset.


## Data sources

Every FRED series has a page at `https://fred.stlouisfed.org/series/<ID>`.

| Variable(s) | Source | Link |
|---|---|---|
| WTI futures, price and volume | Investing.com, saved as `futures1983-2026.csv` | https://www.investing.com/commodities/crude-oil-historical-data |
| WTI spot, monthly | EIA | https://www.eia.gov/dnav/pet/hist/LeafHandler.ashx?n=pet&s=rwtc&f=m |
| WTI spot, daily (DCOILWTICO) | FRED | https://fred.stlouisfed.org/series/DCOILWTICO |
| TB3MS, GS10, BAA, AAA, CPIAUCSL, INDPRO, UNRATE, M2SL, TCU, CFNAI | FRED | https://fred.stlouisfed.org |
| IGREA, Kilian's index of global real economic activity | FRED | https://fred.stlouisfed.org/series/IGREA |
| U.S. crude oil production and imports | EIA (series MCRFPUS2, MCRIMUS2) | https://www.eia.gov/petroleum/data.php |
| U.S. crude oil ending stocks | EIA ([series ID]) | https://www.eia.gov/petroleum/data.php |
| S&P 500, VIX, OVX | Yahoo Finance (^GSPC, ^VIX, ^OVX), downloaded when the script runs | https://finance.yahoo.com |
| EPU | Baker, Bloom and Davis (2016) | https://www.policyuncertainty.com |
| GEPU, MPU, TPU, CPU | Davis (2016); Husted, Rogers and Sun (2020); Caldara et al. (2020); Gavriilidis (2021) | https://www.policyuncertainty.com |
| GPR | Caldara and Iacoviello (2022) | https://www.matteoiacoviello.com/gpr.htm |
| EMV | Baker, Bloom, Davis and Kost (2019) | https://www.policyuncertainty.com |
| IDEMV | Baker, Bloom, Davis, Kost, Sammon and Viratyosin (2020) | https://www.policyuncertainty.com |
| WUI | Ahir, Bloom and Furceri (2022) | https://worlduncertaintyindex.com |
| RA, EU | Bekaert, Engstrom and Xu (2022) | https://www.nancyxu.net/risk-aversion-index |
| RU, MU, FU | Jurado, Ludvigson and Ng (2015); Ludvigson, Ma and Ng (2021) | https://www.sydneyludvigson.com/macro-and-financial-uncertainty-indexes |


## Notes

- S&P 500, VIX and OVX are downloaded each time the script runs, so the
