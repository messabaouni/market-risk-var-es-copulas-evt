# Market Risk Quantification: VaR, ES, Copulas and Extreme Value Theory

Risk quantification project (MESIO UPC-UB, December 2024), written in R.
It measures the market risk of an equally weighted portfolio of Apple (AAPL), Alphabet (GOOG) and Microsoft (MSFT)
using daily closing prices from 9 Sep 2021 to 6 Sep 2024.

## Methods
1. **Value at Risk and Expected Shortfall** with the variance-covariance method (Normal and Student-t) and historical simulation
2. **Multivariate distributions:** normality tests (Jarque-Bera, Shapiro-Wilk, QQ plots) and fitting of Gaussian, t, hyperbolic, NIG, Variance Gamma and generalised hyperbolic models (`ghyp`), compared by AIC, plus Monte Carlo risk estimates
3. **Copulas:** Gumbel and Clayton copulas with Student-t marginals, fitted by maximum pseudo-likelihood and simulated
4. **Extreme Value Theory:** exponential and Pareto family tail fits and GPD / peaks-over-threshold analysis

## Structure
```
data/prices.csv        daily closing prices (AAPL, GOOG, MSFT)
R/risk_analysis.R      full analysis script
report/report.Rmd      report source (includes the script)
report/report.pdf      rendered report
```

## How to read the results
VaR is the loss that is exceeded only with a small probability (for example 1% or 5%) over one day. Expected Shortfall (ES) is the average loss when that threshold is exceeded. The report compares how these estimates change with the model: Normal, Student-t, generalised hyperbolic, copula-based and extreme-value approaches, the latter ones being designed to capture heavy tails. See `report/report.pdf` for the numbers and conclusions.

## Data
Daily closing prices of AAPL, GOOG and MSFT from public market data, as provided for the course.

## How to run
From the repository root:
```r
source("R/risk_analysis.R")
```
Render the report with `rmarkdown::render("report/report.Rmd")`.

Packages: quantmod, fitdistrplus, QRM, pracma, knitr, fBasics, ghyp, copula, MASS, evmix, evir, ercv, EnvStats.

## Author
Mohammed Essabaouni
