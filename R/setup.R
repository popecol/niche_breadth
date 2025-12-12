
# -------------------------------------------------------------------------
# Helper functions & Graphical parameteres
# -------------------------------------------------------------------------

# Accompanies the paper:
# Askarova G., Puchalska E., Lewandowski M., Sexton J., Kuczyński L., Skoracka A. 
# 'Does biotic niche expansion come at an abiotic cost? Evidence from experimental evolution.'


# Libraries ---------------------------------------------------------------

library(mgcv)
library(RColorBrewer)
library(emmeans)


# Graphical parameters ----------------------------------------------------

# display.brewer.all(5, type = "div")
# pal <- rev(brewer.pal(5, "RdBu"))
pal <- brewer.pal(5, "PuOr")
col1 <- pal[1]
col1a <- adjustcolor(col1, alpha.f = 0.01)
col1b <- adjustcolor(col1, alpha.f = 1/4)
col1c <- adjustcolor(col1, alpha.f = 1/2)
col2 <- pal[5]
col2a <- adjustcolor(col2, alpha.f = 0.01)
col2b <- adjustcolor(col2, alpha.f = 1/4)
col2c <- adjustcolor(col2, alpha.f = 1/2)
col <- c(col2, col1)
colb <- c(col2b, col1b)
colc <- c(col2c, col1c)


# Functions ---------------------------------------------------------------

ci <- function(x, alpha = 0.05) {
  # Calculates empirical confidence intervals at the `alpha` level.
  quantile(x, probs = c(alpha / 2, 1 - (alpha / 2)))
}


niche_breadth <- function(x) {
  # Calculates the Levins normalised index of host niche breadth.
  
  # Argument:
  #   x: numeric vector of instantaneous population growth rates (r) across all hosts.
  
  # Returns:
  #   The normalized Levins index value.
  
  # Details:
  #   The output ranges from 0 (no viable growth) to 1 (uniform growth).
  #   If all r-s are <=0, the index is zero.
  
  # The index is 0 if no hosts support positive growth (all r <= 0).
  if(all(x <= 0))
    B <- 0
  
  else {
    # Relative fitness
    p <- x / sum(x)
    
    # The total number of tested hosts
    n <- length(x)
    
    # The index
    B <- 1 / (n * sum(p^2))
  }
  return(B)
}


# Some tests
niche_breadth(c(-10, 1))
niche_breadth(c(0, 1))
niche_breadth(c(0, 10))

niche_breadth(x = c(0, 0, 0, 0, 0, 0))
niche_breadth(x = c(1, 1, 1, 1, 1, 1))
niche_breadth(c(1, 0, 0, 0, 0, 0))
niche_breadth(c(1, 1, 1, 1, 1, 0))

