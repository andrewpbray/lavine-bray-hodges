# # HMO Example
library(dplyr)
library(nlme) # the standard analysis, by REML
library(devtools)
devtools::install_github("andrewpbray/lmmoptim")
library("lmmoptim")

# model 2, with some fixed effects
# random intercept model

data(hmo)
data(hmo_states)
hmobig <- left_join(hmo, hmo_states, by = "state")

y <- hmobig$indPrem
n <- length(y)
mod <- lm(indPrem ~ expPerAdm + (region == "NE") + state, data = hmobig, x = TRUE)
x <- mod$x[, 1:3]
z <- mod$x[, -(1:3), drop = FALSE]
SigE <- diag(n)
SigS <- diag(44)

# the standard reml analysis
mod <- lme(indPrem ~ expPerAdm + (region == "NE"), data = hmobig, random = ~1 | state)
summary(mod)

lines <- findlines(x, z, y, SigE, SigS)

pdf("figs/hmolines_HH11.pdf")
showlines(lines)
dev.off()

mlrebox <- with(lines,
                makebox(lims.sigsqs = c(0, max(int.sigsqs[is.finite(int.sigsqs)])),
	                      lims.sigsqe = c(0, max(int.sigsqe[is.finite(int.sigsqe)])),
	                      status = rep("straddle", nrow(lines)),
                        lines = lines))

boxes.HH11 <- fitlmm(lines = lines, mlrebox, eps = 5, delE = log(10), 
                     delS = log(10), ratio = TRUE, M = 5, maxit = 10)

p <- showboxes(boxes.HH11)
pdf("figs/hmo_HH11_boxes.pdf")
p + scale_x_continuous(name = expression(sigma[e]^2)) +
    scale_y_continuous(name = expression(sigma[s]^2)) +
    theme_bw()
dev.off()

p <- showfunc(boxes.HH11)
pdf("figs/hmo_HH11_rll.pdf")
p + scale_x_log10(name = expression(sigma[e]^2),
                  breaks = c(100,300,1000,3000,10000)) +
    scale_y_log10(name = expression(sigma[s]^2),
                  breaks = c(30,100,300,1000,3000,10000)) +
    theme_bw()
dev.off()

# Run the algorithm again, with input refined as suggested by
# the previous run
box2 <- with(lines,
             makebox(lims.sigsqe = c(300, 1000),
                     lims.sigsqs = c(0, 1000),
	                   status = rep("straddle", nrow(lines)),
                     lines = lines))

boxes.HH11 <- fitlmm(lines = lines, box2, eps = 1, M = 10, maxit = 15)

p <- showboxes(boxes.HH11)
pdf("figs/hmo_HH11_boxes2.pdf")
p + scale_x_continuous(name = expression(sigma[e]^2), limits = c(300, 1000)) +
    scale_y_continuous(name = expression(sigma[s]^2)) +
    theme_bw()
dev.off()

p <- showfunc(boxes.HH11)
pdf("figs/hmo_HH11_rll2.pdf")
p + scale_x_log10(name = expression(sigma[e]^2),
                  breaks = c (300,500,1000)) +
    scale_y_log10(name = expression(sigma[s]^2),
                  breaks = c(1, 3, 10, 30, 100, 300, 1000)) +
    theme_bw()
dev.off()

# Bayesian analysis using Hodges 1998 prior
lines <- addprior(lines, a.E = 1, b.E = 0, a.S = 1.1, b.S = .1)
pdf("figs/hmolines_HH11_Bayes.pdf")
showlines(lines)
dev.off()

boxes.HH11Bayes <- fitlmm(lines = lines, box2, eps = 1, M = 10, maxit = 20)

p <- showboxes(boxes.HH11Bayes)
jpeg("figs/hmo_HH11Bayes_boxes.jpg")
p + scale_x_continuous(name = expression(sigma[e]^2), limits = c(300, 1000)) +
    scale_y_continuous(name = expression(sigma[s]^2)) +
    theme_bw()
dev.off()

p <- showfunc(boxes.HH11Bayes)
jpeg("figs/hmo_HH11Bayes_rll.jpg")
p + scale_x_log10(name = expression(sigma[e]^2),
                  breaks = c(300, 500, 1000)) +
    scale_y_log10(name = expression(sigma[s]^2),
                    # breaks = c (1, 3, 10, 30,100,300,1000)
                  breaks = 10^seq(-3, 3, by = 1)) +
    theme_bw()
dev.off()

# The figure reveals a ridge.  Let's investigate
ridge <- with(boxes.HH11Bayes,
              rll.upper > max(rll.lower) & (rll.upper-rll.lower) > 5)

with(boxes.HH11Bayes, range(sigsqs.lo[good]))
with(boxes.HH11Bayes, range(sigsqs.hi[good]))
with(boxes.HH11Bayes, range(rll.upper[good]))
with(boxes.HH11Bayes, range(rll.lower[good]))
with(boxes.HH11Bayes, range(rll.upper[good] - rll.lower[good]))

tmp <- seq(0, .1, length = 1000)
tmp2 <- dinvgamma(tmp, shape = 1.1, rate = 0.1)
plot(tmp, tmp2, type = "l")
abline(v = 0.04761905)

save.image("model2.RData")

