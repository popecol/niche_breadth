
# -------------------------------------------------------------------------
# Thermal niche breadth
# -------------------------------------------------------------------------

# Accompanies the paper:
# Askarova G., Puchalska E., Lewandowski M., Sexton J., Kuczyński L., Skoracka A. 
# 'Does biotic niche expansion come at an abiotic cost? Evidence from experimental evolution.'


# Setup -------------------------------------------------------------------

source("R/setup.R")


# Data --------------------------------------------------------------------

data <- read.csv("data/thermal_niche.csv", stringsAsFactors = TRUE)

# Calculate the finite growth rate (the initial no. of females was always 5):
data <- transform(data, lambda = n / 5)

summary(data)

lattice::xyplot(lambda ~ temp | regime, data = data, xlab = "Temperature", ylab = "lambda")


# Growth rates at 24 degrees
dat_24 <- subset(data, temp == 24)
lattice::histogram(~ lambda | regime, data = dat_24)
m24 <- lm(lambda ~ regime - 1, dat_24)
summary(m24)
cbind(coef(m24), confint(m24))


# Model ------------------------------------------------------------------

k <- 6
m <- gam(lambda ~ regime + s(temp, by = regime, k = k) + s(temp, rep, bs = "fs", k = k) + s(inv, bs = "re") - 1, family = tw, data = data)

summary(m)
gam.vcomp(m)
anova(m)

# Diagnostics
op <- par(mfrow = c(2, 2)); gam.check(m, rep = 100); par(op)


# Plot the model ----

# The link scale corresponds to instantaneous growth rates:
plot(m, pages = 1, all.terms = TRUE, scale = 0, seWithMean = TRUE, residuals = TRUE, by.resids = TRUE, pch = 21)

# The response scale corresponds to finite growth rates:
plot(m, pages = 1, all.terms = TRUE, scale = 0, seWithMean = TRUE, trans = exp, residuals = TRUE, by.resids = TRUE, pch = 21)


# Simulation from the model -----------------------------------------------

# No. of simulations:
R <- 1e4

# Gaussian approximation:
# br <- rmvn(R, coef(m), vcov(m)) 

# Metropolis Hastings sampler (discard every fifth draw):
br <- gam.mh(m, thin = 5, ns = 5 * R)$bs 
dim(br)

# Predictions from the simulated models
temps <- seq(10, 40, 0.01)
newdat <- expand.grid(regime = levels(data$regime), temp = temps)
Xp <- predict(m, newdat, type = "lpmatrix", exclude = c("s(rep)", "s(inv)", "s(temp,rep)"), newdata.guaranteed = TRUE)
mn <- Xp %*% coef(m)
sim <- Xp %*% t(br)
dim(sim)

idx <- newdat$regime == "Alternating"
sim_alt <- sim[idx, ]
sim_sta <- sim[!idx, ]
mn_alt <- mn[idx]
mn_sta <- mn[!idx]

# Growth rates vs. Temperature ----
# The thermal niche breadth is defined as the range of temperatures at which the instantaneous population growth rate is non-negative.

sidx <- sample(500) # A random sample to speed up drawing.
op <- par(mar = c(4.5, 4.5, 2.5, 2))
matplot(temps, sim_alt[, sidx], type = "l", col = col2a, xlab = "Temperature (\u00B0C)", ylab = "Instantaneous growth rate", cex.lab = 1.3, cex.main = 1.5)
# lines(temps, mn_alt, col = col1)
matlines(temps, sim_sta[, sidx], col = col1a)
# lines(temps, mn_sta, col = col2)
abline(h = 0)
legend("bottomleft", legend = levels(data$regime), lty = 1, col = col, lwd = 2, bty = "n", cex = 1.1)
par(op)

# Optima
opt_alt <- temps[apply(sim_alt, 2, which.max)]
mean(opt_alt)
ci(opt_alt)

opt_sta <- temps[apply(sim_sta, 2, which.max)]
mean(opt_sta)
ci(opt_sta)

# Growth rates at optima
r_alt <- apply(sim_alt, 2, max)
lambda_alt <- exp(r_alt)
mean(lambda_alt); ci(lambda_alt)

r_sta <- apply(sim_sta, 2, max)
lambda_sta <- exp(r_sta)
mean(lambda_sta); ci(lambda_sta)

# Zero crossings (numerical root finding)
zc_alt <- diff(sign(sim_alt))
zc_sta <- diff(sign(sim_sta))

# Lower range
idx_min_alt <- apply(zc_alt, 2, function(x) which(x == 2))
t_min_alt <- (temps[idx_min_alt - 1] + temps[idx_min_alt]) / 2
mean(t_min_alt); ci(t_min_alt)

idx_min_sta <- apply(zc_sta, 2, function(x) which(x == 2))
t_min_sta <- (temps[idx_min_sta - 1] + temps[idx_min_sta]) / 2
mean(t_min_sta); ci(t_min_sta)


# Upper range
idx_max_alt <- apply(zc_alt, 2, function(x) which(x == -2))
t_max_alt <- (temps[idx_max_alt - 1] + temps[idx_max_alt]) / 2
mean(t_max_alt); ci(t_max_alt)

idx_max_sta <- apply(zc_sta, 2, function(x) which(x == -2))
t_max_sta <- (temps[idx_max_sta - 1] + temps[idx_max_sta]) / 2
mean(t_max_sta); ci(t_max_sta)


# Niche breadth
nb_alt <- t_max_alt - t_min_alt
mean(nb_alt); ci(nb_alt)

nb_sta <- t_max_sta - t_min_sta
mean(nb_sta); ci(nb_sta)


# Contrasts ---------------------------------------------------------------

# Differences in optima
d_opt <- opt_alt - opt_sta
hist(d_opt, breaks = 50); ci(d_opt)

# Differences in niche breadth between 'Stable' and 'Alternating'.
delta <- nb_alt - nb_sta
mean(delta); ci(delta)
hist(delta, breaks = 50); ci(delta, 0.01)

