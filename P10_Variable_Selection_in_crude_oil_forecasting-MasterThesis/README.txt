# Variable selection in crude oil forecasting

Code for my master's thesis (MSc Mathematics-Economics, Aalborg University,
2026). The thesis compares penalized regression and boosting methods (LASSO,
Elastic Net, boosting) for selecting predictors when forecasting crude oil
returns and volatility, using recursive out-of-sample forecasts.

The full thesis can be requested for viewing at my personal email.

## Structure

- `Data_creation/` has the R scripts that build the datasets, along with the
  datasets and the raw files they are built from. See its README for the
  sources.
- `ZhangData/` replicates the return forecasting analysis of Zhang et al.
  (2019) on the dataset I built.
- `Volatility/` has the analysis for the volatility dataset.

The thesis contains more analyses than are included here. I left them out to
keep the repo easy to read.
