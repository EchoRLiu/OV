library(nlmixr2)
# library(RxODE)
library(data.table)
library(ggplot2)

model <- function() {
  ini({
    # Fixed effects (initial estimates)
    ttheta1 <- log(0.1)   # rho
    ttheta2 <- log(0.01)  # kappa
    ttheta3 <- log(0.1)   # psi
    ttheta4 <- log(0.1)   # phi
    ttheta5 <- log(0.1)   # alpha
    ttheta6 <- log(0.1)   # beta
    ttheta7 <- log(0.1)   # delta_v

    # Random effects
    eta_Cu0 ~ 0.1         # Variance of random effect on log_Cu0

    # Error model
    prop.err <- 0.1

  })
  model({

    # Transform parameters
    rho <- exp(ttheta1)
    kappa <- exp(ttheta2)
    psi <- exp(ttheta3)
    phi <- exp(ttheta4)
    alpha <- exp(ttheta5)
    beta <- exp(ttheta6)
    delta_v <- exp(ttheta7)

    # Individual initial condition for C_u(0) with random effect
    Cu0 <- exp(log(400) + eta_Cu0)

    # Initial conditions
    C_u(0) <- 1 / (kappa + exp(-rho) * (1/Cu0 - kappa))
    Ci0(0) <- 0
    Ci1(0) <- 0
    Ci2(0) <- 0
    Ci3(0) <- 0
    Ci4(0) <- 0
    Ci5(0) <- 0
    Cl(0)  <- 0
    V(0)   <- 3 * 10^9

    # ODE system
    d/dt(C_u)  = rho * C_u * (1 - kappa * (C_u + Ci0 + Ci1 + Ci2 + Ci3 + Ci4 + Ci5 + Cl)) - psi * V * C_u
    d/dt(Ci0)  = rho * Ci0 * (1 - kappa * (C_u + Ci0 + Ci1 + Ci2 + Ci3 + Ci4 + Ci5 + Cl)) + psi * V * C_u - phi * Ci0
    d/dt(Ci1)  = rho * Ci1 * (1 - kappa * (C_u + Ci0 + Ci1 + Ci2 + Ci3 + Ci4 + Ci5 + Cl)) + phi * Ci0 - phi * Ci1
    d/dt(Ci2)  = rho * Ci2 * (1 - kappa * (C_u + Ci0 + Ci1 + Ci2 + Ci3 + Ci4 + Ci5 + Cl)) + phi * Ci1 - phi * Ci2
    d/dt(Ci3)  = rho * Ci3 * (1 - kappa * (C_u + Ci0 + Ci1 + Ci2 + Ci3 + Ci4 + Ci5 + Cl)) + phi * Ci2 - phi * Ci3
    d/dt(Ci4)  = rho * Ci4 * (1 - kappa * (C_u + Ci0 + Ci1 + Ci2 + Ci3 + Ci4 + Ci5 + Cl)) + phi * Ci3 - phi * Ci4
    d/dt(Ci5)  = rho * Ci5 * (1 - kappa * (C_u + Ci0 + Ci1 + Ci2 + Ci3 + Ci4 + Ci5 + Cl)) + phi * Ci4 - phi * Ci5
    d/dt(Cl)   = rho * Cl  * (1 - kappa * (C_u + Ci0 + Ci1 + Ci2 + Ci3 + Ci4 + Ci5 + Cl)) + phi * Ci5 - alpha * Cl
    d/dt(V)    = beta * alpha * Cl - psi * V * C_u - delta_v * V

    # Observation model: Observable is TotalCells
    TotalObs <- C_u + Ci0 + Ci1 + Ci2 + Ci3 + Ci4 + Ci5 + Cl
    TotalObs ~ prop(prop.err)
  })
}

# Load the data
data <- read.table("measurements.tsv", header=TRUE, sep="\t")

# Extract the relevant columns and create the data frame
ID <- seq_len(nrow(data))  # Create an id column as a unique identifier for each row
TIME <- data$time   # Use the 'time' column for TIME
DV <- data$measurement  # Use 'measurement' column for DV
Group <- data$simulationConditionId  # Use 'simulationConditionId' for Group

# Create the final data frame
data_final <- data.frame(ID=ID, TIME=TIME, DV=DV, Group=Group)

# Display the first few rows
head(data_final)

# Convert data to data frame
data_final <- as.data.frame(data_final)

# Convert Group to factor (if not already)
data_final$Group <- as.factor(data_final$Group)

head(data_final)

# Fit the model
fit <- nlmixr2(model, data=data_final, est="saem", list(print=0), table=list(cwres=TRUE, npde=TRUE))

# Summarize the fit
print(summary(fit))

# Plot the fit
plot(fit)