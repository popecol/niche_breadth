
# -------------------------------------------------------------------------
# Figure 2. Host niche breadth
# -------------------------------------------------------------------------

# Accompanies the paper:
# Askarova G., Puchalska E., Lewandowski M., Sexton J., Kuczyński L., Skoracka A. 
# 'Does biotic niche expansion come at an abiotic cost? Evidence from experimental evolution.'

# Note: Run 'biotic_niche.R' first.


# Prepare the data for the graph --------------------------------------------

# Estimated population means and their standard errors
pre <- predict(m, newdat, exclude = c("s(rep)", "s(inv)"), se = TRUE, newdata.guaranteed = TRUE)
fit <- pre$fit; se <- pre$se.fit

# Asymptotic confidence intervals
lci <- fit + 1.96 * se; uci <- fit - 1.96 * se
p <- data.frame(fit, lci, uci)
# p[] <- exp(p) # Uncomment to get finite growth rates.
newd <- data.frame(newdat, p)

# Expected niche breadth across regimes
niche_breadth(fit[seq(1, 12, 2)]) # Alternating
niche_breadth(fit[seq(2, 12, 2)]) # Stable

# Cooking the x-axis values:
nhosts <- nlevels(data$host)
remove_every <- 3
ng <- nhosts * remove_every
at1 <- (1:ng)[-c(remove_every * 1:nhosts)]
at2 <- seq(remove_every / 2, ng, remove_every)

# Labels
labels <- levels(newd$regime)

levels(newd$host)
hosts <- c("Oat-grass", "Brome", "Quackgrass", "Barley", "Rye", "Wheat")

# Kernel density estimates
d_nb_sta <- density(nb_sta, bw = 0.01); plot(d_nb_sta)
d_nb_alt <- density(nb_alt, bw = d_nb_sta$bw); plot(d_nb_alt)
dx <- mean(nb_sta)
d <- density(dx + delta, bw = d_nb_sta$bw, from = from)

# The range of the X-axis (set manually to improve eligibility)
from <- -0.04; to <- 0.65


# The Graph ---------------------------------------------------------------

mult <- 1

# Dimensions in millimetres:
w <- 173; h <- 173 * 2/3

# Dimensions in inches * multiplier:
wi <- round(mult * w / 25.4, 1)
hi <- round(mult * h / 25.4, 1)

cex.lab <- 1.7
cex.axis <- 1.5

op <- par(no.readonly = TRUE)

cairo_pdf("Figures/Fig_2.pdf", width = wi, height = hi, symbolfamily = "OpenSymbol")

nf <- layout(matrix(c(1, 2, 2, 1, 3, 3), 2, 3, byrow = TRUE))
# layout.show(nf)


# Fig. 2A -----------------------------------------------------------------
# The instantaneous population growth rates of experimental populations evolved in 'Stable' and 'Alternating' host environments on the six plant species that were tested. Vertical bars denote 95% confidence intervals.

par(mar = c(10, 5, 3, 1))

plot(at1, newd$fit, type = "h", lwd = 6, lend = 1, xaxt = "n", xlab = "", ylab = "Growth rate", cex.lab = 1.9, xlim = c(0.5, 17.5), ylim = range(p), col = colc)
abline(h = 0, col = "grey30")
segments(at1, newd$lci, at1, newd$uci, col = col, lwd = 1)
axis(1, at = at2, labels = hosts, cex.axis = cex.axis, las = 2)
# mtext("Host plant", 1, 4, cex = 2)
legend(2.3, -1.5, legend = rev(labels), fill = rev(colc), cex = 1.5, bty = "n", text.width = 4, ncol = 1)
mtext("(a)", side = 2, line = -1, padj = -7.7, cex = 1.5, las = 1)


# Fig. 2B ----
# The posterior distributions of the host niche breadth index estimated for the 'Stable' and 'Alternating' regimes. The vertical lines represent the estimated means.

par(mar = c(5, 6, 3, 1))

plot(d_nb_sta, xlab = "Host niche breadth index", ylab = "", cex.lab = cex.lab, main = "", xlim = c(from, to), ylim = range(d_nb_sta$y, d_nb_alt$y))
polygon(d_nb_sta, col = col1b, border = NA)
abline(v = mean(nb_sta))
lines(d_nb_alt)
polygon(d_nb_alt, col = col2b, border = NA)
abline(v = mean(nb_alt))
legend(0.15, 12, legend = rev(labels), fill = rev(colc), cex = 1.5, bty = "n", text.width = 4, ncol = 1)
mtext("(b)", side = 2, line = -1, padj = -3.7, cex = 1.5, las = 1)

mtext("Probability density", side = 2, line = 3, adj = 2.5, cex = 1.3)


# Fig. 2C ----
# The posterior distributions of the effect size (Δ), which is calculated as the difference between the estimated niche breadths for the 'Alternating' and 'Stable' regimes. The thin vertical line represents the mean effect size and the dashed and dotted lines represent the 95% and 99% empirical confidence intervals, respectively. The thick vertical line at zero is a reference: if niche breadths were not significantly different at the given alpha level, the confidence intervals for Δ would include zero.

par(mar = c(5, 6, 2, 1))

plot(d, xlab = expression(paste("Effect size (", Delta, ")")), ylab = "", cex.lab = cex.lab, main = "", xlim = c(from, to), xaxt = "n")
polygon(d, col = adjustcolor("grey", alpha.f = 1/4), border = NA)
abline(v = dx, lwd = 2)
abline(v = c(mean(dx + delta), ci(dx + delta)), lty = c(1, 2, 2))
abline(v = ci(dx + delta, 0.01), lty = 3)
lab <- seq(-0.1, 1, 0.1)
axis(1, las = 1, at = lab + dx, labels = lab)
mtext("(c)", side = 2, line = -1, padj = -4, cex = 1.5, las = 1)

par(op)

dev.off()
