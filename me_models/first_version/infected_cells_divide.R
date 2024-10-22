library(nlmixr2)
library(rxode2)
library(data.table)
library(ggplot2)
library(reshape2)

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
    add.err <- 0.1
    
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
    Cu0 <- exp(log(400)+ eta_Cu0)
    
    # Initial conditions
    C_u(0) <- 1 / (kappa + exp(-rho) * (1/Cu0 - kappa))
    Ci0(0) <- 0
    Ci1(0) <- 0
    Ci2(0) <- 0
    Ci3(0) <- 0
    Ci4(0) <- 0
    Ci5(0) <- 0
    Cl(0)  <- 0
    
    # Use the preprocessed V0_INIT column from the dataset
    V(0) <- V0_INIT
    
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
    logTotalObs <- log(TotalObs)
    
    logTotalObs ~ add(add.err)
  })
}

# Load the data from the CSV file
df <- read.csv("/Users/yuhongliu/Documents/OV/data/ov_datasets_v7.csv")

# Increase the number of digits after the dot for data.frames
options(digits=10)
# Extract tumor volume data for vvDD condition
tumor_vol_vvdd <- df[2:6, seq(2, 20, by = 2)]  # Select every other column starting from the second column (Y values)
# Create a DataFrame with known time points and extracted tumor volumes
tumor_vol_vvdd_df <- data.frame(lapply(tumor_vol_vvdd, as.numeric))
colnames(tumor_vol_vvdd_df) <- paste0("vvdd_sample_", 1:10)
rownames(tumor_vol_vvdd_df) <- 0:4

# Extract tumor volume data for normal/PBS condition
tumor_vol_pbs <- df[2:6, seq(22, 40, by = 2)]  # Select every other column starting from the 22nd column (Y values)
# Create a DataFrame with known time points and extracted tumor volumes
tumor_vol_pbs_df <- data.frame(lapply(tumor_vol_pbs, as.numeric))
colnames(tumor_vol_pbs_df) <- paste0("pbs_sample_", 1:10)
rownames(tumor_vol_pbs_df) <- 0:4

# Convert data to numeric type
tumor_vol_vvdd_df[] <- lapply(tumor_vol_vvdd_df, as.numeric)
tumor_vol_pbs_df[] <- lapply(tumor_vol_pbs_df, as.numeric)

# Define time points
time_points <- 0:4

# Add time points to both dataframes
tumor_vol_vvdd_df$TIME <- time_points
tumor_vol_pbs_df$TIME <- time_points

# Reshape vvDD dataframe
vvdd_melted <- melt(tumor_vol_vvdd_df, id.vars = "TIME", 
                    variable.name = "ID", value.name = "DV")
vvdd_melted$ID <- as.numeric(gsub("vvdd_sample_", "", vvdd_melted$ID))  # Extract sample number from ID
vvdd_melted$GROUP <- "vvDD"  # Add GROUP column for vvDD

# Reshape PBS dataframe
pbs_melted <- melt(tumor_vol_pbs_df, id.vars = "TIME", 
                   variable.name = "ID", value.name = "DV")
pbs_melted$ID <- as.numeric(gsub("pbs_sample_", "", pbs_melted$ID)) + 10  # Sample number + 10
pbs_melted$GROUP <- "ctrl"  # Add GROUP column for PBS

# Combine both dataframes
data_final <- rbind(vvdd_melted, pbs_melted)

# Display the first few rows
head(data_final)

# Assign initial values of V(0) based on GROUP
data_final$V0_INIT <- ifelse(data_final$GROUP == "vvDD", 3 * 10^9, 0)

# Log-transform the observed tumor volume data (DV)
data_final$DV <- log(data_final$DV + 1)  # Add a small constant to handle zeros if applicable

# Convert data to data frame
data_final <- as.data.frame(data_final)

# Convert Group to factor (if not already)
data_final$GROUP <- as.factor(data_final$GROUP)

head(data_final)

# Fit the model
fit <- nlmixr2(model, data=data_final, est="saem", list(print=0), table=list(cwres=TRUE, npde=TRUE))

# Summarize the fit
print(summary(fit))

# Extract fitted values
fit_data <- as.data.frame(fit)

# Plot DV against PRED for each five rows in one plot
plot_list <- list()
for (i in seq(1, nrow(fit_data), by = 5)) {
  plot_data <- fit_data[i:(i+4), ]
  p <- ggplot(plot_data, aes(x = TIME)) +
    geom_point(aes(y = DV), color = "blue") +
    geom_line(aes(y = PRED), color = "red") +
    geom_point(aes(y = PRED), color = "red") +
    labs(title = paste("Rows", i, "to", i+4), x = "Time", y = "Value") +
    theme_minimal()
  plot_list[[length(plot_list) + 1]] <- p
}

# Display the plots one by one
for (i in seq_along(plot_list)) {
  print(plot_list[[i]])
  Sys.sleep(1)  # Pause for 1 second between plots
}

# Save each plot to a separate file
# for (i in seq_along(plot_list)) {
#  ggsave(filename = paste0("plot_", i, ".png"), plot = plot_list[[i]], width = 8, height = 6)
#}