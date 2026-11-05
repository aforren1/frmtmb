# What brms 2.23.0 builds for an nl body under the families that have no
# single mu: categorical, multinomial, dirichlet, dirichlet_multinomial,
# logistic_normal, and mixtures. stancode only. Seed 20261102.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
set.seed(20261102)
n <- 100
d <- data.frame(x = rnorm(n), z = rnorm(n), ntr = 10L)
d$cat <- factor(sample(c("a", "b", "c"), n, TRUE))
d$cnt <- t(rmultinom(n, 10, c(1, 1, 1)))
colnames(d$cnt) <- c("a", "b", "c")
d$sim <- d$cnt / 10 * 0.8 + 0.2 / 3
d$yg <- rnorm(n)
show <- function(lab, f, fam, data = d) {
  r <- tryCatch(suppressWarnings(suppressMessages(
    stancode(f, data = data, family = fam))), error = function(e) e)
  cat("\n==", lab, "\n")
  if (inherits(r, "error")) {
    cat("  REFUSES:", gsub("\n", " ", conditionMessage(r)), "\n")
    return(invisible())
  }
  ln <- strsplit(r, "\n")[[1]]
  keep <- grep(paste0("nlp_|mu[a-zA-Z0-9_]*\\[n\\] =|mu[a-zA-Z0-9_]* =|",
                      "mu[a-zA-Z0-9_]* \\+=|target \\+= [a-z_]+_l"), ln,
               value = TRUE)
  cat(paste0("  ", unique(trimws(keep))), sep = "\n")
}
show("categorical nl=TRUE", bf(cat ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
     categorical())
show("categorical nlf(mub), mub-only lf",
     bf(cat ~ 1, nlf(mub ~ b0 + b1 * x), b0 ~ 1, b1 ~ 1, nl = TRUE),
     categorical())
show("categorical mub ~ nlpars, flag",
     bf(cat ~ 1) + nlf(mub ~ b0 + b1 * x) + lf(b0 ~ 1) + lf(b1 ~ 1) +
       lf(muc ~ x), categorical())
show("multinomial nl=TRUE",
     bf(cnt | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
     multinomial())
show("dirichlet nl=TRUE", bf(sim ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
     dirichlet())
show("dirichlet_multinomial nl=TRUE",
     bf(cnt | trials(ntr) ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
     dirichlet_multinomial())
show("logistic_normal nl=TRUE", bf(sim ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
     logistic_normal())
show("mixture nl=TRUE", bf(yg ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE),
     mixture(gaussian(), gaussian()))
show("mixture nlf(mu1)",
     bf(yg ~ 1) + nlf(mu1 ~ b0 + b1 * x) + lf(b0 ~ 1) + lf(b1 ~ 1) +
       lf(mu2 ~ x), mixture(gaussian(), gaussian()))
show("mvbf nl", bf(yg ~ b0 + b1 * x, b0 + b1 ~ 1, nl = TRUE) + bf(z ~ x) +
       set_rescor(FALSE), gaussian())
