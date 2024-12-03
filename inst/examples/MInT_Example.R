# Set the seed for reproducibility
set.seed(123)

# Create a sequence of 120 numbers
x <- seq(1, 120)

# Generate the columns
AA <- sin(x*pi/6) + rnorm(120, 0, 1) # Sine component with random noise
AB <- 0.05*x + rnorm(120, 0, 0.5) # Linear component
B <- cos(x*pi/6)+ rnorm(120, 0, 1) # Cosine component

# Combine the columns into a matrix
matrix <- cbind(AA, AB, B)
hts = ts(matrix, frequency = 12)

# Define S matrix
S <- rbind(c(1,1,1), c(1,1,0), diag(1,3))
rownames(S) <-c("Total", "A", "AA", "AB", "B")
colnames(S) <- c("AA", "AB", "B")

# Aggregate hts on all levels
hts.complete <- ts(t(S %*% t(hts)), frequency = 12)

# Fit a model to the time series
hts.models = lapply(hts.complete, function(c.ts) forecast::ets(c.ts))

# Fit a model to the time series
hts.models = lapply(hts.complete, function(c.ts) forecast::ets(c.ts))
# Generate predictions based on this model
hts.forecasts = sapply(hts.models, function(mdl) forecast::forecast(mdl, h = 1)$mean)
# Extract residuals
hts.residuals = sapply(hts.models, function(mdl) mdl$residuals)

# Compute reconciled forecasts
MinT(fcasts = hts.forecasts, Smat = S, residual = hts.residuals)
