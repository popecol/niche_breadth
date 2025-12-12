
# -------------------------------------------------------------------------
# Host niche breadth 
# -------------------------------------------------------------------------

# Accompanies the paper:
# Askarova G., Puchalska E., Lewandowski M., Sexton J., Kuczyński L., Skoracka A. 
# 'Does biotic niche expansion come at an abiotic cost? Evidence from experimental evolution.'


# Setup -------------------------------------------------------------------

source("R/setup.R")


# Data --------------------------------------------------------------------

data <- read.csv("data/host_niche.csv", stringsAsFactors = TRUE)

# Calculate the finite growth rate (the initial no. of females was always 5):
data <- transform(data, lambda = n / 5)

summary(data)

lattice::bwplot(lambda ~ host | regime, data = data, xlab = "Host", ylab = "Growth rate")


# Model ------------------------------------------------------------------

# Cell means model:
m <- gam(lambda ~ regime:host + s(rep, bs = "re") + s(inv, bs = "re") - 1, family = tw, data = data)

summary(m)
gam.vcomp(m)
anova(m)

# Diagnostics
op <- par(mfrow = c(2, 2)); gam.check(m, rep = 100); par(op)
plot.gam(m, pages = 1, all.terms = TRUE)

emmeans(m, pairwise ~ regime | host, type = "response")
emmip(m, ~ regime | host, type = "response", CIs = TRUE)

# The model could be re-parameterised to test the main factors and interactions separately.
# ma <- gam(lambda ~ regime * host + s(rep, bs = "re") + s(inv, bs = "re") - 1, family = tw, data = data)
# summary(ma)
# anova(ma)


# Simulation from the model -----------------------------------------------

# No. of simulations:
R <- 1e4

# Gaussian approximation:
# br <- rmvn(R, coef(m), vcov(m)) 

# Metropolis Hastings sampler (discard every fifth draw):
br <- gam.mh(m, thin = 5, ns = 5 * R)$bs 
dim(br)

# Predictions from simulated models
newdat <- expand.grid(regime = levels(data$regime), host = levels(data$host))
Xp <- predict(m, newdat, type = "lpmatrix", exclude = c("s(rep)", "s(inv)"), newdata.guaranteed = TRUE)
mn <- Xp %*% coef(m)
sim <- Xp %*% t(br)
dim(sim)

idx <- newdat$regime == "Alternating"
sim_alt <- sim[idx, ]
sim_sta <- sim[!idx, ]
mn_alt <- mn[idx]
mn_sta <- mn[!idx]

nb_alt <- apply(sim_alt, 2, niche_breadth)
nb_sta <- apply(sim_sta, 2, niche_breadth)
xlim <- range(nb_alt, nb_sta)

op <- par(mfrow = c(2, 1))
hist(nb_alt, xlim = xlim, breaks = 50, main = "Alternating regime", xlab = "Host niche breadth")
abline(v = c(mean(nb_alt), ci(nb_alt)), lty = c(1, 2, 2), lwd = c(2, 1, 1))
hist(nb_sta, xlim = xlim, breaks = 50, main = "Stable regime", xlab = "Host niche breadth")
abline(v = c(mean(nb_sta), ci(nb_sta)), lty = c(1, 2, 2), lwd = c(2, 1, 1))
par(op)


# Contrasts ---------------------------------------------------------------

# The effect size (Δ) was calculated as the difference between the niche breadth estimated for 'Alternating' and 'Stable' regimes. Thus, it shows how much the niche breadth of 'Alternating' regime differs from the niche breadth of 'Stable' regime.

# To test whether the effect size was significantly different from zero, we derived the distribution of Δ by simulating (10,000 of times) from a posterior distribution of the GAMM parameters.

# Then, empirical 95% confidence intervals for Δ were calculated by finding 0.025 and 0.975 quantiles of the sampling distribution.

# Differences in niche breadth between 'alternated' and 'stable'.
delta <- nb_alt - nb_sta
mean(delta)
hist(delta, breaks = 50); ci(delta)

# To reproduce Figure 2, go to the 'Fig_2.R' script.
