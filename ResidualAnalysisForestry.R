## Residual analysis on Forestry dataset using area as the target variable

# Read data
forestry <- read.csv("forestry.csv", header = T)

# Full model
model <- lm(area ~ height + caliper + htcal, data = forestry)
summ <- summary(model) # Full model summary

# Derive Studentized residuals
n <- dim(forestry)[1]
sigma_hat <- summ$sigma # Sigma-hat
X <- cbind(rep(1, n), forestry$height, forestry$caliper, forestry$htcal) # X matrix
H <- X %*% solve(t(X) %*% X) %*% t(X) # Hat matrix
h <- diag(H) # Diagonal elements of H
st_resid <- model$residuals/(sigma_hat*sqrt(1-h)) # Studentized Residuals

par(mfrow = c(2,2))

# Plot of Studentized Residuals vs Index
index <- 1:n
{plot(st_resid, ylim = c(min(-3, min(st_resid)),max(3, max(st_resid))), 
     col = "darkblue", xlab = "Index", pch = 16,
     ylab = "Studentized Residuals", main = "Studentized Residuals vs Index")
abline(h=c(0,3,-3), lty = c(1,2,2), col = "red")}

# Plot of Studentized Residuals vs Fitted Values
{plot(model$fitted.values, st_resid, ylim = c(min(-3, min(st_resid)), max(3, max(st_resid))), 
     col = "darkblue", pch = 16, xlab = "Fitted Values", ylab = "Studentized Residuals", 
     main = "Studentized Residuals vs Fitted Values")
abline(h=c(0,3,-3), lty = c(1,2,2), col = "red")}

# Histogram of Studentized Residuals
hist(x = st_resid, xlab = "Studentized Residuals", 
     main = "Histogram of Studentized Residuals", breaks=seq(-5,5,0.25))

# QQ-Plot of Studentized Residuals
{qqnorm(st_resid, pch = 16, main = "QQ-plot of Studentized Residuals")
qqline((st_resid), col = "red", lwd = 2)}

# Get observation with largest studentized residual
larStRes <- max(st_resid)
indexLarStRes <- which(st_resid == larStRes)

# Leverage calculations
lev <- hatvalues(model)
n <- length(lev)
index <- 1:n
avLev = mean(lev)
highLev = which(lev > (2*avLev)) # Observations with leverage considered to be high

# Cook's D-statistic calculations
cooks <- cooks.distance(model)
n <- length(cooks)
index <- 1:n
sortedCooks <- sort(cooks, decreasing = T)
highestCooks <- sortedCooks[1:3]
mostInf <- c()
for(i in 1:3) { # Get indices of top 3 most influential observations
  mostInf <- append(mostInf, which(cooks == highestCooks[i]))
}

par(mfrow = c(1,2))

# Plot of Leverage vs Index
{plot(x = index, y = lev, col = "purple", pch = 16, xlab = "Index", 
     ylab = "Hat Values", ylim = c(0,1), main = "Leverage vs Index")
abline(h = (2*avLev), col = "red", lty = 2)}

# Plot of Cook's D-statistic vs Index
{plot(x = index, y = cooks, col = "blue", pch = 16, xlab = "Index", 
     ylab = "Cook's D-statistic", ylim = c(0,1), main = "Cook's D-statistic vs Index")
abline(h = 0.5, col = "red", lty = 2)}

# Delete observations 10 and 29 since they look like outliers
newForestry <- forestry[-c(10,29),]

# New model
newModel <- lm(area ~ height + caliper + htcal, data = newForestry)
newSumm <- summary(newModel) # Summary of new model

# Derive Studentized residuals
n <- dim(newForestry)[1]
sigma_hat <- newSumm$sigma # Sigma-hat
X <- cbind(rep(1, n), newForestry$height, newForestry$caliper, newForestry$htcal) # X matrix
H <- X %*% solve(t(X) %*% X) %*% t(X) # Hat matrix
h <- diag(H) # Diagonal elements of H
stRes <- newModel$residuals/(sigma_hat*sqrt(1-h)) # Studentized Residuals

par(mfrow = c(2,2))

# Plot of Studentized Residuals vs Index
index <- 1:n
{plot(stRes, ylim = c(min(-3, min(stRes)),max(3, max(stRes))),col = "purple", pch = 16, 
      xlab = "Index", ylab = "Studentized Residuals", main = "Studentized Residuals vs Index")
abline(h=c(0,3,-3), lty = c(1,2,2), col = "red")}

# Plot of Studentized Residuals vs Fitted Values
{plot(newModel$fitted.values, stRes, ylim = c(min(-3, min(stRes)),max(3, max(stRes))),
     col = "darkblue", pch = 16, 
     xlab = "Fitted Values", ylab = "Studentized Residuals", 
     main = "Studentized Residuals vs Fitted Values")
abline(h=c(0,3,-3), lty = c(1,2,2), col = "red")}

# Histogram of Studentized Residuals
hist(x = stRes, xlab = "Studentized Residuals", 
     main = "Histogram of Studentized Residuals", breaks=seq(-5,5,0.25))

# QQ-Plot of Studentized Residuals
{qqnorm(stRes, pch = 16, main = "QQ-plot of Studentized Residuals")
qqline((stRes), col = "red", lwd = 2)}

## Perform log transformation to see if variability stabilises
# Model of log(area) without observarions 10 and 29
logModel <- lm(log(area) ~ height + caliper + htcal, data = newForestry)
logSumm <- summary(logModel)

# Derive Studentized residuals
n <- dim(newForestry)[1]
sigma_hat <- logSumm$sigma # Sigma-hat
X <- cbind(rep(1, n), newForestry$height, newForestry$caliper, newForestry$htcal) # X matrix
H <- X %*% solve(t(X) %*% X) %*% t(X) # Hat matrix
h <- diag(H) # Diagonal elements of H
logStRes <- logModel$residuals/(sigma_hat*sqrt(1-h)) # Studentized Residuals

# Plot of Studentized Residuals vs Fitted Values
{plot(logModel$fitted.values, logStRes, ylim = c(min(-3, min(logStRes)),max(3, max(logStRes))), 
     col = "darkblue", pch = 16, xlab = "Fitted Values", ylab = "Studentized Residuals", 
     main = "Log-transformation")
abline(h=c(0,3,-3), lty = c(1,2,2), col = "red")}


