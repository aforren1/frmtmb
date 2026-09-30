# Reviewer, claim 2: gr(g, cov = A) on a subsetted response whose rows
# lack a level of A, against the separate fit on its rows, and brms's
# level count. Seed 224. Log: dev/aterms2-rev-log-02d-cov.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(e) suppressWarnings(suppressMessages(e))
r <- function(e) tryCatch(q(e), error = function(e) paste("ERROR:", conditionMessage(e)))
set.seed(224)
n <- 180
G <- 9
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(letters[1:G], length.out = n)))
A <- 0.5 ^ abs(outer(1:G, 1:G, "-"))
dimnames(A) <- list(letters[1:G], letters[1:G])
u <- as.vector(t(chol(A)) %*% rnorm(G))
d$y1 <- 1 + d$x + u[d$g] + rnorm(n)
d$y2 <- d$z + rnorm(n)
d$s <- d$g != "i"
fj <- r(frm(bf(y1 | subset(s) ~ x + (1 | gr(g, cov = A))) + bf(y2 ~ z),
            data = d, data2 = list(A = A), family = gaussian()))
fs <- r(frm(y1 ~ x + (1 | gr(g, cov = A)), data = d[d$s, ], data2 = list(A = A)))
f2 <- q(frm(y2 ~ z, data = d))
cat("joint:", if (is.character(fj)) fj else format(logLik(fj), digits = 13), "\n")
cat("separate y1 on its rows:", if (is.character(fs)) fs else
  format(as.numeric(logLik(fs)) + as.numeric(logLik(f2)), digits = 13), "\n")
sd <- tryCatch(brms::standata(brms::bf(y1 | subset(s) ~ x + (1 | gr(g, cov = A))) +
                                brms::bf(y2 ~ z) + brms::set_rescor(FALSE),
                              data = d, data2 = list(A = A)),
               error = function(e) conditionMessage(e))
cat("brms N_1:", if (is.list(sd)) sd$N_1 else sd, "\n")
# does a level of A with no rows get marginalized (A[-i, -i]) or
# conditioned on zero? The first equals brms, which keeps the level
A8 <- A[1:8, 1:8]
fs8 <- r(frm(y1 ~ x + (1 | gr(g, cov = A8)), data = droplevels(d[d$s, ]),
             data2 = list(A8 = A8)))
cat("separate with A[-i,-i]:", format(logLik(fs8), digits = 13),
    " with full A:", format(logLik(fs), digits = 13), "\n")
