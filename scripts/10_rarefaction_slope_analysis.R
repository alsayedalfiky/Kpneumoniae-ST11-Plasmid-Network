# Rarefaction Slope and Saturation Analysis
# Quantifying the rate of phylogenetic discovery and curve saturation

# 1. Calculate Discovery Rates (Slopes) at the end of the curves
# Note: Step sizes are 50 for Total Population and 20 for ST11 as per Step 09
final_idx_total <- nrow(df_total)
slope_total <- (df_total$Mean_PD[final_idx_total] - df_total$Mean_PD[final_idx_total - 1]) / 50

final_idx_st11 <- nrow(df_st11)
slope_st11 <- (df_st11$Mean_PD[final_idx_st11] - df_st11$Mean_PD[final_idx_st11 - 1]) / 20

# 2. Calculate Saturation Coefficients (Final Slope / Initial Slope)
# A lower coefficient indicates the curve has reached a plateau
initial_slope_total <- (df_total$Mean_PD[2] - df_total$Mean_PD[1]) / 50
sat_total <- slope_total / initial_slope_total

initial_slope_st11 <- (df_st11$Mean_PD[2] - df_st11$Mean_PD[1]) / 20
sat_st11 <- slope_st11 / initial_slope_st11

# 3. Output Analysis Summary
cat("\n--- Saturation and Slope Metrics ---\n")
cat("Total Population - Final Discovery Rate:", round(slope_total, 6), "\n")
cat("ST11 Lineage     - Final Discovery Rate:", round(slope_st11, 6), "\n")
cat("Total Population - Saturation Coefficient:", round(sat_total, 4), "\n")
cat("ST11 Lineage     - Saturation Coefficient:", round(sat_st11, 4), "\n")
