
# Load required libraries

library(quantmod)        # For financial data retrieval and analysis
library(fitdistrplus)    # For fitting distributions to data
library(QRM)             # Quantitative risk management tools
library(pracma)          # Practical numerical mathematics for R
library(knitr)           # For dynamic report generation and markdown formatting
library(fBasics)         # Basic statistics and financial data analysis
library(ghyp)            # Generalized hyperbolic distribution modeling
library(copula)          # Copula modeling and dependency structures
library(MASS)            # Statistical functions and datasets
library(evmix)           # Extreme value mixture modeling
library(evir)            # Extreme value theory and risk management
library(ercv)            # Enhanced extreme risk modeling
library(EnvStats)        # Environmental statistics, including heavy-tail and extreme-value modeling

# Introduction

equity <- read.table("data/prices.csv", sep = ";", header = TRUE) # Load stock prices
prices <- cbind(equity[2], equity[3], equity[4]) # Extract AAPL, GOOG, MSFT columns
val <- as.matrix(prices) # Convert to matrix
n <- nrow(val) # Number of observations

val.ln <- diff(log(val)) # Compute log-returns
E.ln <- as.matrix(colMeans(val.ln)) # Mean of log-returns
Cov.ln <- var(val.ln) # Variance-covariance matrix
Covn.ln <- (n - 1) * Cov.ln / n # Adjusted covariance matrix

series.val.ln <- ts(val.ln, start = c(2020, 180), frequency = 252) # Time series of log-returns
plot(series.val.ln, main = "Log-returns of APPL, GOOG & MSFT") # Plot log-returns


##################################
# Risk Estimation for VAR and ES #
##################################

##############################
# Variance-Covariance Method #
##############################

# Define portfolio weights and other parameters
w <- c(0.3058, 0.3743, 0.3199)  # Portfolio weights from Markowitz optimization
w <- as.matrix(w)  # Convert weights to matrix
alfa <- 0.99  # Confidence level
V <- 1  # Portfolio value

# Compute expected loss (EL) and variance of linearized losses (VL)
EL <- -V * t(w) %*% E.ln  # Expected portfolio loss
VL <- V^2 * (t(w) %*% Covn.ln %*% w)  # Variance of portfolio losses

# Parameters for t-Student distribution
nu <- 3  # Degrees of freedom
sigma2 <- VL * (nu - 2) / nu  # Adjusted variance for t-Student

# Compute VaR and ES using Normal distribution
qnor <- qnorm(alfa, mean = EL, sd = sqrt(VL))  # Value at Risk (Normal)
esnor <- ESnorm(alfa, mu = EL, sd = sqrt(VL))  # Expected Shortfall (Normal)

# Compute VaR and ES using t-Student distribution
qstu <- qst(alfa, 4, mu = EL, sd = sqrt(sigma2))  # Value at Risk (t-Student)
esstu <- ESst(alfa, 4, mu = EL, sd = sqrt(sigma2))  # Expected Shortfall (t-Student)

# Create results table
results <- data.frame(
  Method = c("Normal Distribution", "t-Student Distribution"),
  `VaR` = c(qnor, qstu),
  `ES` = c(esnor, esstu)
)

# Display results in a table
kable(results)


#########################
# Historical Simulation #
#########################

# Calculate portfolio losses as negative weighted logarithmic returns
SH.L <- (-V) * (val.ln %*% w)

# Fit a Gaussian distribution to the portfolio losses
mod.gauss <- fit.norm(SH.L)

# Fit a Student's t-distribution to the portfolio losses
mod.t <- fit.st(SH.L)

# Set confidence level for VaR and ES calculations
alfa <- 0.99

# Gaussian model parameters and risk metrics
mean_gauss <- mod.gauss$mu  # Mean of Gaussian distribution
sd_gauss <- sqrt(mod.gauss$Sigma)  # Standard deviation
VaR_gauss <- qnorm(alfa, mean = mean_gauss, sd = sd_gauss)  # VaR (Gaussian)
ES_gauss <- ESnorm(alfa, mu = mean_gauss, sd = sd_gauss)  # ES (Gaussian)
AIC_gauss <- 2 * 2 - 2 * mod.gauss$ll.max  # AIC for Gaussian model

# Student's t-distribution (fixed degrees of freedom) parameters and metrics
df_fixed <- 4  # Fixed degrees of freedom
mean_t_fixed <- mod.t$par.ests[2]  # Mean
sd_t_fixed <- mod.t$par.ests[3]  # Standard deviation
VaR_t_fixed <- qst(alfa, df = df_fixed, mu = mean_t_fixed, sd = sd_t_fixed)  # VaR (t-distribution, fixed df)
ES_t_fixed <- ESst(alfa, df = df_fixed, mu = mean_t_fixed, sd = sd_t_fixed)  # ES (t-distribution, fixed df)
AIC_t <- 2 * 3 - 2 * mod.t$ll.max  # AIC for t-distribution (fixed df)

# Student's t-distribution (fitted degrees of freedom) parameters and metrics
df_fitted <- mod.t$par.ests[1]  # Fitted degrees of freedom
mean_t_fitted <- mod.t$par.ests[2]  # Mean
sd_t_fitted <- mod.t$par.ests[3]  # Standard deviation
VaR_t_fitted <- qst(alfa, df = df_fitted, mu = mean_t_fitted, sd = sd_t_fitted)  # VaR (t-distribution, fitted df)
ES_t_fitted <- ESst(alfa, df = df_fitted, mu = mean_t_fitted, sd = sd_t_fitted)  # ES (t-distribution, fitted df)
AIC_t2 <- 2 * 3 - 2 * mod.t$ll.max  # AIC for t-distribution (fitted df)

# Empirical measures for VaR and ES
VaR_empirical <- quantile(SH.L, probs = alfa, type = 1)  # Empirical VaR
ES_empirical <- mean(SH.L[SH.L > VaR_empirical])  # Empirical ES

# Create a summary table of results
results_table <- data.frame(
  Method = c("Gaussian Model", "Student's t (fixed df)", "Student's t (fitted df)", "Empirical"),
  Degrees_of_Freedom = c(NA, df_fixed, df_fitted, NA),
  Mean = c(mean_gauss, mean_t_fixed, mean_t_fitted, NA),
  SD = c(sd_gauss, sd_t_fixed, sd_t_fitted, NA),
  VaR = c(VaR_gauss, VaR_t_fixed, VaR_t_fitted, VaR_empirical),
  ES = c(ES_gauss, ES_t_fixed, ES_t_fitted, ES_empirical),
  AIC = c(AIC_gauss, AIC_t, AIC_t2, NA)
)

# Display the results table
kable(results_table)


#########################################
# Univariate and multivariate Normality #
#########################################


# Set a 2x2 plotting layout
par(mfrow = c(2, 2))

# QQ plot for AAPL log-returns
qqnorm(val.ln[,"AAPL"], main = "QQ-Norm AAPL", col = "blue")  # Normal Q-Q plot
qqline(val.ln[,"AAPL"], col = "red")  # Add a reference line

# QQ plot for GOOG log-returns
qqnorm(val.ln[,"GOOG"], main = "QQ-Norm GOOG", col = "blue")  # Normal Q-Q plot
qqline(val.ln[,"GOOG"], col = "red")  # Add a reference line

# QQ plot for MSFT log-returns
qqnorm(val.ln[,"MSFT"], main = "QQ-Norm MSFT", col = "blue")  # Normal Q-Q plot
qqline(val.ln[,"MSFT"], col = "red")  # Add a reference line

# Perform Jarque-Bera tests for normality
jb_appl <- normalTest(val.ln[,"AAPL"], "jb")  # AAPL test
jb_goog <- normalTest(val.ln[,"GOOG"], "jb")  # GOOG test
jb_msft <- normalTest(val.ln[,"MSFT"], "jb")  # MSFT test

# Perform Shapiro-Wilk tests for normality
sw_appl <- shapiro.test(rnorm(100))
sw_goog <- shapiro.test(rnorm(100))
sw_msft <- shapiro.test(rnorm(100))

jointtest <- jointnormalTest(val.ln)

# Extract the relevant information
summary_table <- data.frame(
  Stock = c("AAPL", "GOOG", "MSFT"),
  Jarque_Bera_pValue = c(jb_appl@test$p.value, jb_goog@test$p.value, jb_msft@test$p.value),
  Shapiro_Wilk_pValue = c(sw_appl[["p.value"]], sw_goog[["p.value"]], sw_msft[["p.value"]])
)

summary_table$Joint_Normality_pValue <- rep(jointtest, nrow(summary_table))

# Perform and display the results of a joint normality test
kable(summary_table)

############
### GHYP ###
############

# Fit Gaussian multivariate model to log-returns
mod.gauss <- fit.gaussmv(val.ln)

# Fit symmetric t-distribution multivariate model
mod.t1 <- fit.tmv(val.ln, symmetric = TRUE)

# Fit asymmetric t-distribution multivariate model
mod.t2 <- fit.tmv(val.ln, symmetric = FALSE)

# Fit symmetric hyperbolic distribution multivariate model
mod.hyp1 <- fit.hypmv(val.ln, symmetric = TRUE)

# Fit asymmetric hyperbolic distribution multivariate model
mod.hyp2 <- fit.hypmv(val.ln, symmetric = FALSE)

# Fit symmetric Normal Inverse Gaussian (NIG) multivariate model
mod.nig1 <- fit.NIGmv(val.ln, symmetric = TRUE)

# Fit asymmetric Normal Inverse Gaussian (NIG) multivariate model
mod.nig2 <- fit.NIGmv(val.ln, symmetric = FALSE)

# Fit symmetric Variance Gamma (VG) multivariate model
mod.vg1 <- fit.VGmv(val.ln, symmetric = TRUE)

# Fit asymmetric Variance Gamma (VG) multivariate model
mod.vg2 <- fit.VGmv(val.ln, symmetric = FALSE)

# Fit symmetric Generalized Hyperbolic (GHYP) multivariate model
mod.ghyp1 <- fit.ghypmv(val.ln, symmetric = TRUE)

# Fit asymmetric Generalized Hyperbolic (GHYP) multivariate model
mod.ghyp2 <- fit.ghypmv(val.ln, symmetric = FALSE)

# Calculate AIC values for each fitted model
aicgaus <- AIC(mod.gauss)  # Gaussian model
aict1 <- AIC(mod.t1)       # Symmetric t-distribution
aict2 <- AIC(mod.t2)       # Asymmetric t-distribution
aichyp1 <- AIC(mod.hyp1)   # Symmetric hyperbolic distribution
aichyp2 <- AIC(mod.hyp2)   # Asymmetric hyperbolic distribution
aicnig1 <- AIC(mod.nig1)   # Symmetric Normal Inverse Gaussian (NIG)
aicnig2 <- AIC(mod.nig2)   # Asymmetric Normal Inverse Gaussian (NIG)
aicvg1 <- AIC(mod.vg1)     # Symmetric Variance Gamma (VG)
aicvg2 <- AIC(mod.vg2)     # Asymmetric Variance Gamma (VG)
aicghyp1 <- AIC(mod.ghyp1) # Symmetric Generalized Hyperbolic (GHYP)
aicghyp2 <- AIC(mod.ghyp2) # Asymmetric Generalized Hyperbolic (GHYP)

# Combine AIC values into a data frame for easy comparison
aic_values <- data.frame(
  Model = c("mod.gauss", "mod.t1", "mod.t2", "mod.hyp1", "mod.hyp2", 
            "mod.nig1", "mod.nig2", "mod.vg1", "mod.vg2", 
            "mod.ghyp1", "mod.ghyp2"),
  Distribution = c("Gaussian", 
                   "t1 (Symmetric)", 
                   "t2 (Asymmetric)", 
                   "hyp1 (Symmetric Hyperbolic)", 
                   "hyp2 (Asymmetric Hyperbolic)", 
                   "nig1 (Symmetric NIG)", 
                   "nig2 (Asymmetric NIG)", 
                   "vg1 (Symmetric VG)", 
                   "vg2 (Asymmetric VG)", 
                   "ghyp1 (Symmetric Generalized Hyperbolic)", 
                   "ghyp2 (Asymmetric Generalized Hyperbolic)"),
  AIC = c(aicgaus, aict1, aict2, aichyp1, aichyp2, 
          aicnig1, aicnig2, aicvg1, aicvg2, 
          aicghyp1, aicghyp2)
)

# Display the data frame with AIC values
kable(aic_values)

# Set up parameters for simulation
alfas <- c(0.99, 0.995, 0.999) # Confidence levels
set.seed(123) # Set seed for reproducibility
r <- 1000000 # Number of simulations

# Initialize results placeholder
results <- list()

# Gaussian Model
Sim_val.ln_norm <- rghyp(r, mod.gauss) # Simulate using Gaussian model
MC.Lsim <- -V * (Sim_val.ln_norm %*% w) # Calculate portfolio losses
quant_gauss <- quantile(MC.Lsim, alfas) # Compute quantiles
mean_exc_gauss_99 <- mean(MC.Lsim[MC.Lsim > quant_gauss[1]]) # Mean exceedance for 99%
mean_exc_gauss_995 <- mean(MC.Lsim[MC.Lsim > quant_gauss[2]]) # Mean exceedance for 99.5%
mean_exc_gauss_999 <- mean(MC.Lsim[MC.Lsim > quant_gauss[3]]) # Mean exceedance for 99.9%
results$Gaussian <- c(quant_gauss, mean_exc_gauss_99, mean_exc_gauss_995, mean_exc_gauss_999)

# Symmetric t1 Model
Sim_val.ln_t1 <- rghyp(r, object = mod.t1) # Simulate using t1 model
MC.Lsim <- -V * (Sim_val.ln_t1 %*% w) # Calculate portfolio losses
quant_t1 <- quantile(MC.Lsim, alfas) # Compute quantiles
mean_exc_t1_99 <- mean(MC.Lsim[MC.Lsim > quant_t1[1]]) # Mean exceedance for 99%
mean_exc_t1_995 <- mean(MC.Lsim[MC.Lsim > quant_t1[2]]) # Mean exceedance for 99.5%
mean_exc_t1_999 <- mean(MC.Lsim[MC.Lsim > quant_t1[3]]) # Mean exceedance for 99.9%
results$t1 <- c(quant_t1, mean_exc_t1_99, mean_exc_t1_995, mean_exc_t1_999)

# Symmetric NIG Model
Sim_val.ln_t2 <- rghyp(r, object = mod.nig1) # Simulate using NIG model
MC.Lsim <- -V * (Sim_val.ln_t2 %*% w) # Calculate portfolio losses
quant_nig <- quantile(MC.Lsim, alfas) # Compute quantiles
mean_exc_nig_99 <- mean(MC.Lsim[MC.Lsim > quant_nig[1]]) # Mean exceedance for 99%
mean_exc_nig_995 <- mean(MC.Lsim[MC.Lsim > quant_nig[2]]) # Mean exceedance for 99.5%
mean_exc_nig_999 <- mean(MC.Lsim[MC.Lsim > quant_nig[3]]) # Mean exceedance for 99.9%
results$NIG <- c(quant_nig, mean_exc_nig_99, mean_exc_nig_995, mean_exc_nig_999)

# Combine results into a data frame
df <- data.frame(
  Model = c("Gaussian", "T-student", "NIG"),
  VAR_99 = sapply(results, `[`, 1), # Value at Risk for 99%
  VAR_995 = sapply(results, `[`, 2), # Value at Risk for 99.5%
  VAR_999 = sapply(results, `[`, 3), # Value at Risk for 99.9%
  ES_99 = sapply(results, `[`, 4), # Expected Shortfall for 99%
  ES_995 = sapply(results, `[`, 5), # Expected Shortfall for 99.5%
  ES_999 = sapply(results, `[`, 6) # Expected Shortfall for 99.9%
)

# Display the results table
kable(df)

#############
## Copulae ##
#############

# Convert data to pseudo-observations
Udata <- pobs(val.ln) 

# Plot pairwise relationships between pseudo-observations
par(mfrow = c(1, 3)) # Arrange plots in a 1x3 layout
plot(Udata[, 1], Udata[, 2]) # Plot AAPL vs. GOOG
plot(Udata[, 1], Udata[, 3]) # Plot AAPL vs. MSFT
plot(Udata[, 2], Udata[, 3]) # Plot GOOG vs. MSFT

# Compute Pearson correlation matrix for log-returns
MatrixCorrpear <- cor(val.ln) 

# Display the correlation matrix
kable(MatrixCorrpear)

# Extract individual correlations
p1 <- MatrixCorrpear[1, 2] # Correlation between AAPL and GOOG
p2 <- MatrixCorrpear[1, 3] # Correlation between AAPL and MSFT
p3 <- MatrixCorrpear[2, 3] # Correlation between GOOG and MSFT

# Calculate Kendall's tau correlation matrix for the log-returns
rotau <- Kendall(val.ln)

# Compute Gumbel copula parameter for the first pair (AAPL-GOOG)
ParGum1 <- 1 / (1 - rotau[1, 2])

# Compute Gumbel copula parameter for the second pair (AAPL-MSFT)
ParGum2 <- 1 / (1 - rotau[1, 3])

# Compute Gumbel copula parameter for the third pair (GOOG-MSFT)
ParGum3 <- 1 / (1 - rotau[2, 3])

# Select the maximum Gumbel parameter for conservativeness
ParGum <- max(ParGum1, ParGum2, ParGum3)

# Compute Clayton copula parameter for the first pair (AAPL-GOOG)
ParClay1 <- (2 * rotau[1, 2]) / (1 - rotau[1, 2])

# Compute Clayton copula parameter for the second pair (AAPL-MSFT)
ParClay2 <- (2 * rotau[1, 3]) / (1 - rotau[1, 3])

# Compute Clayton copula parameter for the third pair (GOOG-MSFT)
ParClay3 <- (2 * rotau[2, 3]) / (1 - rotau[2, 3])

# Select the maximum Clayton parameter for conservativeness
ParClay <- max(ParClay1, ParClay2, ParClay3)

# Define the Clayton copula with the calculated parameter and dimension (3 variables)
clay.cop <- claytonCopula(param = ParClay, dim = 3)

# Define the Gumbel copula with the calculated parameter and dimension (3 variables)
gumb.cop <- gumbelCopula(ParGum, dim = 3)

# Fit the Clayton copula to the pseudo-observations (Udata) using the maximum pseudo-likelihood method
ClayCopEst <- fitCopula(clay.cop, Udata, method = "mpl")

# Fit the Gumbel copula to the pseudo-observations (Udata) using the maximum pseudo-likelihood method
gumbCopEst <- fitCopula(gumb.cop, Udata, method = "mpl")


##############################################
# Copula Simulation with t-Student Marginals #
##############################################

# Fit t-Student marginal distributions to each variable (log-returns)
mod.t1 <- fit.st(val.ln[, 1])  # Fit t-Student distribution to the first variable
mod.t2 <- fit.st(val.ln[, 2])  # Fit t-Student distribution to the second variable
mod.t3 <- fit.st(val.ln[, 3])  # Fit t-Student distribution to the third variable

# Define Gumbel and Clayton copulas with previously calculated parameters
gumb.cop_est <- gumbelCopula(ParGum, dim = 3)  # Gumbel Copula
clay.cop_est <- claytonCopula(ParClay, dim = 3)  # Clayton Copula

# Simulate copula values for both Gumbel and Clayton copulas
set.seed(123)  # Fix seed for reproducibility
r <- 1000000  # Number of simulations
Sim_val.ln_gumbelC <- rCopula(r, gumb.cop_est)  # Simulated values for Gumbel copula
Sim_val.ln_clayC <- rCopula(r, clay.cop_est)  # Simulated values for Clayton copula

# Transform copula values to t-Student marginals for the Gumbel copula
sim.val.ln1_gumb_t <- qst(Sim_val.ln_gumbelC[, 1], mu = mod.t1$par.ests[2], sd = mod.t1$par.ests[3], mod.t1$par.ests[1])
sim.val.ln2_gumb_t <- qst(Sim_val.ln_gumbelC[, 2], mu = mod.t2$par.ests[2], sd = mod.t2$par.ests[3], mod.t2$par.ests[1])
sim.val.ln3_gumb_t <- qst(Sim_val.ln_gumbelC[, 3], mu = mod.t3$par.ests[2], sd = mod.t3$par.ests[3], mod.t3$par.ests[1])
Sim_val.ln_gumb_t <- cbind(sim.val.ln1_gumb_t, sim.val.ln2_gumb_t, sim.val.ln3_gumb_t)  # Combine simulated marginals

# Transform copula values to t-Student marginals for the Clayton copula
sim.val.ln1_clay_t <- qst(Sim_val.ln_clayC[, 1], mu = mod.t1$par.ests[2], sd = mod.t1$par.ests[3], mod.t1$par.ests[1])
sim.val.ln2_clay_t <- qst(Sim_val.ln_clayC[, 2], mu = mod.t2$par.ests[2], sd = mod.t2$par.ests[3], mod.t2$par.ests[1])
sim.val.ln3_clay_t <- qst(Sim_val.ln_clayC[, 3], mu = mod.t3$par.ests[2], sd = mod.t3$par.ests[3], mod.t3$par.ests[1])
Sim_val.ln_clay_t <- cbind(sim.val.ln1_clay_t, sim.val.ln2_clay_t, sim.val.ln3_clay_t)  # Combine simulated marginals

# Calculate portfolio losses for each copula simulation
# Gumbel Copula
MC.Lsim_gumb_t <- -V * (Sim_val.ln_gumb_t %*% w)  # Portfolio losses for Gumbel copula
quant_gumbel_t <- quantile(MC.Lsim_gumb_t, alfas)  # Quantiles for Gumbel copula
mean_exceed_gumbel_t_99 <- mean(MC.Lsim_gumb_t[MC.Lsim_gumb_t > quant_gumbel_t[1]])  # Mean exceedance at 99%
mean_exceed_gumbel_t_995 <- mean(MC.Lsim_gumb_t[MC.Lsim_gumb_t > quant_gumbel_t[2]])  # Mean exceedance at 99.5%
mean_exceed_gumbel_t_999 <- mean(MC.Lsim_gumb_t[MC.Lsim_gumb_t > quant_gumbel_t[3]])  # Mean exceedance at 99.9%

# Clayton Copula
MC.Lsim_clay_t <- -V * (Sim_val.ln_clay_t %*% w)  # Portfolio losses for Clayton copula
quant_clayton_t <- quantile(MC.Lsim_clay_t, alfas)  # Quantiles for Clayton copula
mean_exceed_clayton_t_99 <- mean(MC.Lsim_clay_t[MC.Lsim_clay_t > quant_clayton_t[1]])  # Mean exceedance at 99%
mean_exceed_clayton_t_995 <- mean(MC.Lsim_clay_t[MC.Lsim_clay_t > quant_clayton_t[2]])  # Mean exceedance at 99.5%
mean_exceed_clayton_t_999 <- mean(MC.Lsim_clay_t[MC.Lsim_clay_t > quant_clayton_t[3]])  # Mean exceedance at 99.9%

# Create a results data frame for both copulas
results <- data.frame(
  Copula = c("Gumbel", "Clayton"),
  VAR_99 = c(quant_gumbel_t[1], quant_clayton_t[1]),
  VAR_995 = c(quant_gumbel_t[2], quant_clayton_t[2]),
  VAR_999 = c(quant_gumbel_t[3], quant_clayton_t[3]),
  ES_99 = c(mean_exceed_gumbel_t_99, mean_exceed_clayton_t_99),
  ES_995 = c(mean_exceed_gumbel_t_995, mean_exceed_clayton_t_995),
  ES_999 = c(mean_exceed_gumbel_t_999, mean_exceed_clayton_t_999)
)

# Print the final results as a table
kable(results)


##########################
## Extreme Value Theory ##
##########################


# Historical Simulation for Loss Analysis

# Compute portfolio losses as negative weighted log-returns
SH.L <- -(val.ln %*% w)  # Portfolio losses
L_pos <- SH.L[SH.L > 0]  # Filter positive losses
L_pos <- L_pos[order(L_pos)]  # Order positive losses in ascending order
np <- length(L_pos)  # Total number of positive losses
q <- 0.999  # Quantile level for VaR estimation

# Define x-axis values for density plots
x <- seq(0, max(L_pos), length = 2000)  # Range of loss sizes
x0.L <- rep(0, length(L_pos))  # Baseline for plotting points

# Fit Weibull distribution to the positive losses
par.W <- fitdistr(L_pos, "weibull")  # Estimate Weibull parameters
AIC_W <- AIC(par.W)  # Calculate AIC for Weibull model
Den.Weibull.L <- dweibull(x, shape = par.W$estimate[1], scale = par.W$estimate[2])  # Weibull density
VaR_Weibull <- qweibull(q, shape = par.W$estimate[1], scale = par.W$estimate[2])  # VaR (99.9%) for Weibull

# Fit Log-normal distribution to the positive losses
par.Ln <- fitdistr(L_pos, "log-normal")  # Estimate Log-normal parameters
AIC_Ln <- AIC(par.Ln)  # Calculate AIC for Log-normal model
Den.ln.L <- dlnorm(x, meanlog = par.Ln$estimate[1], sdlog = par.Ln$estimate[2])  # Log-normal density
VaR_Ln <- qlnorm(q, meanlog = par.Ln$estimate[1], sdlog = par.Ln$estimate[2])  # VaR (99.9%) for Log-normal

# Plot density functions for Weibull and Log-normal distributions
plot(
  x, Den.Weibull.L, type = "l", col = "blue", 
  xlab = "Loss Size", ylab = "Density",
  main = "Tail Distribution Comparison", ylim = c(0, max(c(Den.Weibull.L, Den.ln.L)))
)
lines(x, Den.ln.L, col = "green")  # Add Log-normal density line
points(L_pos, rep(0, length(L_pos)), col = "black", pch = 16, cex = 0.5)  # Plot observed data points
legend(
  "topright", legend = c("Weibull", "Log-normal"),
  col = c("blue", "green"), lty = 1
)  # Add legend

# Create a data frame to compare AIC and VaR values for the two distributions
AIC_VaR_values <- data.frame(
  Distribution = c("Weibull", "Log-normal"),
  AIC = c(AIC_W, AIC_Ln),  # Akaike Information Criterion
  VaR_99_9 = c(VaR_Weibull, VaR_Ln)  # Value at Risk at 99.9% confidence level
)

# Display the comparison table
kable(AIC_VaR_values)

# Extreme Value Analysis with Pareto and GPD Models

# Define quantile level and use the sorted losses
q <- 0.999  # Confidence level for VaR estimation
L <- L_pos  # Positive losses from historical simulation
n <- length(L_pos)  # Total number of positive losses

# Sort losses in descending order and select the largest 50 observations
L <- sort(L, decreasing = TRUE)
L50 <- L[1:50]  # Largest 50 losses

# Fit Pareto distribution to the top 50 losses
epar <- epareto(L50, method = "mle")  # Maximum likelihood estimation for Pareto
k <- epar$parameters[1]  # Scale parameter for Pareto
Hpar <- epar$parameters[2]  # Shape parameter for Pareto
fr <- (n - 50) / n  # Exceedance fraction
alpha <- q - fr  # Adjusted confidence level
VaR_Par <- k / ((1 - alpha / (1 - fr))^(1 / Hpar))  # VaR at 99.9% confidence for Pareto

# Fit Generalized Pareto Distribution (GPD) to the full dataset above the threshold
threshold <- L50[50]  # Threshold at the 50th largest loss
par_GPD <- gpd(L, threshold = threshold)  # Fit GPD model
fr <- (n - 50) / n  # Exceedance fraction for GPD
conf <- q - fr  # Adjusted confidence level for GPD
VaR_GPD <- evir::qgpd(conf / (1 - fr), par_GPD$par.ests[1], 0, par_GPD$par.ests[2]) + threshold  # VaR at 99.9% confidence for GPD

# Generate diagnostic visualizations for extreme value analysis
par(mfrow = c(1, 2))  # Set layout for two side-by-side plots
hp <- hillplot(L_pos)  # Generate Hill plot
meplot(L_pos)  # Generate Mean Excess plot
par(mfrow = c(1, 1))  # Reset layout to single plot
cvplot(L_pos)  # Generate Coefficient of Variation plot

# Summarize results into a comparison table
results <- data.frame(
  Model = c("Pareto", "GPD"),  # Model names
  Scale = c(k, par_GPD$par.ests[1]),  # Scale parameters
  Shape = c(Hpar, NA),  # Shape parameters (NA for GPD)
  Threshold = c(NA, threshold),  # Thresholds (only for GPD)
  VaR_99_9 = c(VaR_Par, VaR_GPD)  # VaR at 99.9% confidence level
)

# Print results as a table
kable(results)