## F2 follow-up: validate_prior() row order against brms's, with priors
## each model takes, and the full rows of base vs lane (saved per arm).
## Seed 21 (data as in dev/review-brmsnames-fr1.R).
##   Rscript dev/review-brmsnames-fr1vp.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
set.seed(21)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), xp = runif(n, 0.5, 2),
                k = factor(sample(1:3, n, TRUE)), k0 = sample(0:3, n, TRUE),
                g = factor(rep(1:12, 20)), t = rep(1:20, each = 12),
                f = factor(sample(c("a", "b", "c"), n, TRUE)))
d$y <- 1 + d$x + 0.3 * d$z^2 + 0.5 * log(d$xp) + rnorm(12)[d$g] + rnorm(n)
d$y2 <- d$x + rnorm(n)
d$yo <- cut(d$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
models <- list(
  list(f = "y ~ exp(x) + z + (1 + z | g)", fam = "gaussian",
       p = list(list("normal(0, 2)", class = "b", coef = "expx"),
                list("exponential(1)", class = "sd", group = "g"))),
  list(f = "bf(y ~ x + f, sigma ~ z)", fam = "gaussian",
       p = list(list("normal(0, 2)", class = "b", dpar = "sigma"),
                list("normal(0, 1)", class = "b", coef = "fb"))),
  list(f = paste("mvbf(bf(y ~ x + (1 | g)), bf(y2 ~ log(xp) + (1 | g)),",
                 "rescor = FALSE)"), fam = "gaussian",
       p = list(list("normal(0, 2)", class = "b", resp = "y2"),
                list("exponential(1)", class = "sd", resp = "y"))),
  list(f = "bf(y ~ a + b * x, a ~ exp(z) + (1 | g), b ~ 1 + z, nl = TRUE)",
       fam = "gaussian",
       p = list(list("normal(0, 2)", class = "b", nlpar = "b"),
                list("normal(1, 2)", class = "b", nlpar = "a"))),
  list(f = "bf(yo ~ x + cs(z))", fam = "acat",
       p = list(list("normal(0, 2)", class = "b"))),
  list(f = "y ~ x + ar(time = t, gr = g)", fam = "gaussian",
       p = list(list("normal(0, 2)", class = "b"),
                list("normal(0, 0.5)", class = "ar")))
)
key <- function(p) with(p, paste(class, coef, group, resp, dpar, nlpar,
                                 sep = "/"))
out <- list()
for (m in models) {
  fr <- eval(parse(text = m$f))
  fb <- eval(parse(text = gsub("(mvbf|bf)[(]", "brms::\\1(", m$f)))
  fam <- get(m$fam, asNamespace("frmtmb"))()
  pf <- do.call(c, lapply(m$p, function(a) do.call(set_prior, a)))
  pb <- Reduce(`+`, lapply(m$p, function(a) do.call(brms::set_prior, a)))
  vp <- tryCatch(as.data.frame(q(validate_prior(pf, fr, data = d,
                                                family = fam))),
                 error = function(e) conditionMessage(e))
  bv <- tryCatch(as.data.frame(q(brms::validate_prior(pb, fb, data = d,
                                                      family = m$fam))),
                 error = function(e) conditionMessage(e))
  dp <- as.data.frame(q(default_prior(fr, data = d, family = fam)))
  out[[m$f]] <- list(dp = dp, vp = vp)
  res <- if (is.character(vp)) paste("frmtmb ERROR", substr(vp, 1, 100))
    else if (is.character(bv)) paste("brms ERROR", substr(bv, 1, 100))
    else {
      ka <- key(vp); kb <- key(bv)
      identical(ka[ka %in% kb], kb[kb %in% ka])
    }
  cat(sprintf("%-70s validate_prior order = brms: %s\n", m$f, res))
}
saveRDS(out, sprintf("dev/review-brmsnames-out/fr1vp-tables-%s.rds", arm))
