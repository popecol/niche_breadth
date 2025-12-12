
# -------------------------------------------------------------------------
# Figure 3. Thermal niche breadth
# -------------------------------------------------------------------------

# Accompanies the paper:
# Askarova G., Puchalska E., Lewandowski M., Sexton J., Kuczyński L., Skoracka A. 
# 'Does biotic niche expansion come at an abiotic cost? Evidence from experimental evolution.'

# Note: Run 'thermal_niche.R' first.


# Prepare the data for the graph --------------------------------------------

nd <- expand.grid(regime = levels(data$regime), temp = seq(10, 40, len = 100))
pre <- predict(m, nd, exclude = c("s(rep)", "s(inv)", "s(temp,rep)"), se = TRUE, newdata.guaranteed = TRUE)
fit <- pre$fit; se <- pre$se.fit

# Asymptotic confidence intervals
lci <- fit + 1.96 * se; uci <- fit - 1.96 * se
p <- data.frame(fit, lci, uci)
# p[] <- exp(p) # Uncomment to get finite growth rates.
nd <- data.frame(nd, p)

labels <- levels(data$regime)

ylim <- range(nd$lci, nd$uci)

# Kernel density estimates
from <- floor(min(nb_alt, nb_sta) * 2) / 2
# to <- ceiling(max(nb_alt, nb_sta) * 2) / 2
to <- 26.5

d_nb_alt <- density(nb_alt, from = from, to = to, bw = 0.1)
d_nb_sta <- density(nb_sta, from = from, to = to, bw = d_nb_alt$bw)

dx <- mean(nb_sta)
d <- density(dx + delta, bw = d_nb_sta$bw, from = from)


# Graphs ------------------------------------------------------------------

mult <- 1

# Dimensions in millimetres:
w <- 173; h <- 173 * 2/3

# Dimensions in inches * multiplier:
wi <- round(mult * w / 25.4, 1)
hi <- round(mult * h / 25.4, 1)

cex.lab <- 1.7
cex.axis <- 1.5

op <- par(no.readonly = TRUE)

cairo_pdf("Figures/Fig_3.pdf", width = wi, height = hi, symbolfamily = "OpenSymbol")

nf <- layout(matrix(c(1, 2, 2, 1, 3, 3), 2, 3, byrow = TRUE))
# layout.show(nf)


# Fig. 3A -----------------------------------------------------------------

par(mar = c(5.2, 5, 3, 1))

plot(fit ~ temp, nd, xlim = c(10, 40), ylim = ylim, type = "n", xlab = "Temperature (\u00B0C)", ylab = "Growth rate", cex.lab = 1.9, main = "")
abline(h = 0, col = "grey30")

for (i in seq(2)) {
  ndd <- subset(nd, regime == labels[i])
  polygon(c(ndd$temp, rev(ndd$temp)), c(ndd$lci, rev(ndd$uci)), col = colb[i], border = NA)
 lines(fit ~ temp, ndd, col = col[i])
}

legend(9, -7, legend = rev(labels), fill = rev(colc), cex = 1.5, bty = "n", text.width = 4, ncol = 1)
mtext("(a)", side = 2, line = -1, padj = -9.3, cex = 1.5, las = 1)


# Fig. 3B -----------------------------------------------------------------

par(mar = c(5, 6, 3, 1))

plot(d_nb_sta, xlab = "Thermal niche breadth (\u00B0C)", ylab = "", cex.lab = cex.lab, main = "", xlim = c(from, to), ylim = range(d_nb_sta$y, d_nb_alt$y))
polygon(d_nb_sta, col = col1b, border = NA)
abline(v = mean(nb_sta))
lines(d_nb_alt)
polygon(d_nb_alt, col = col2b, border = NA)
abline(v = mean(nb_alt))
legend("topright", legend = labels, fill = c(col2b, col1b), cex = 1.5, inset = c(-0.06, -0.4), xpd = TRUE, horiz = TRUE, bty = "n", text.width = 1.7)
mtext("(b)", side = 2, line = -1, padj = -3.7, cex = 1.5, las = 1)
mtext("Probability density", side = 2, line = 3, adj = 2.5, cex = 1.3)


# Fig. 3C -----------------------------------------------------------------

par(mar = c(5, 6, 2, 1))

plot(d, xlab = expression(paste("Effect size (", Delta, ")")), ylab = "", cex.lab = cex.lab, main = "", xlim = c(from, to), xaxt = "n")
polygon(d, col = adjustcolor("grey", alpha.f = 1/4), border = NA)
abline(v = dx, lwd = 2)
abline(v = c(mean(dx + delta), ci(dx + delta, 0.05)), lty = c(1, 2, 2))
abline(v = ci(dx + delta, 0.01), lty = 3)
# summary(delta)
lab <- seq(-4, 1, 1)
axis(1, las = 1, at = lab + dx, labels = lab)
mtext("(c)", side = 2, line = -1, padj = -4, cex = 1.5, las = 1)

par(op)

dev.off()

