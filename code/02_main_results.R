# Load dataset
library(haven)
library(data.table)
library(ivreg)
library(dplyr)

data <- read_dta("C:/Users/eugen/OneDrive - UC San Diego/Desktop/UC San Diego/2024/Fall 2024/Computation 280/Replication/data/data_reg_R.dta")

# Regression models for psi using fixest (feols)
model1 <- feols(infl_reg_sign ~ L4_mean_une + L4_rp | constant, data = data, cluster = "statecode")
model2 <- feols(infl_reg_sign ~ L4_mean_une + L4_rp | statecode, data = data, cluster = "statecode")
model3 <- feols(infl_reg_sign ~ L4_mean_une + L4_rp | statecode + date, data = data, cluster = "statecode")

# Instrumental Variable Regression for psi using ivreg
iv_model <- ivreg(infl_reg_sign ~ L4_mean_une + L4_rp | L4_d20_qt_bartik_sa + L4_rp, data = data)