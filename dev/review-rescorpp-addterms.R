## Review rescorpp: addition terms and autocor under rescor; does the
## joint draw drop what the per-response simulator honored?
##   Rscript dev/review-rescorpp-addterms.R lane
source("C:/Users/adf44/source/r/frmtmb-wt-rescorpp/dev/review-rescorpp-common.R")
suppressMessages(loadNamespace("brms"))
set.seed(71); n <- 120
x <- rnorm(n)
E <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, .7, .7, 1), 2))
d <- data.frame(x = x, y1 = 1 + x + E[, 1], y2 = E[, 2], w = runif(n, 0.5, 2),
                tt = seq_len(n), s1 = 0.5)
d$y1t <- pmax(d$y1, 0.01)
d$ce <- ifelse(d$y1 > 2, "right", "none")
cases <- list(
  trunc = list(b = function(B) B$mvbf(B$bf(y1t | trunc(lb = 0) ~ x), B$bf(y2 ~ x)) + B$set_rescor(TRUE)),
  cens = list(b = function(B) B$mvbf(B$bf(y1 | cens(ce) ~ x), B$bf(y2 ~ x)) + B$set_rescor(TRUE)),
  weights = list(b = function(B) B$mvbf(B$bf(y1 | weights(w) ~ x), B$bf(y2 | weights(w) ~ x)) + B$set_rescor(TRUE)),
  se = list(b = function(B) B$mvbf(B$bf(y1 | se(s1, sigma = TRUE) ~ x), B$bf(y2 ~ x)) + B$set_rescor(TRUE)),
  ar = list(b = function(B) B$mvbf(B$bf(y1 ~ x + ar(tt)), B$bf(y2 ~ x + ar(tt))) + B$set_rescor(TRUE)),
  arcov = list(b = function(B) B$mvbf(B$bf(y1 ~ x + ar(tt, cov = TRUE)), B$bf(y2 ~ x + ar(tt, cov = TRUE))) + B$set_rescor(TRUE)),
  mo = list(b = function(B) B$mvbf(B$bf(y1 ~ x), B$bf(y2 ~ x, sigma ~ x)) + B$set_rescor(TRUE)))
Bb <- asNamespace("brms")
Bf <- asNamespace("frmtmb")
for (nm in names(cases)) {
  cat("\n####", nm, "\n")
  sb <- tryCatch({brms::make_stancode(cases[[nm]]$b(Bb), data = d); "brms stancode OK"},
                 error = function(e) paste("brms ERROR:", conditionMessage(e)))
  cat(sb, "\n")
  fit <- tryCatch(quiet(frm(cases[[nm]]$b(Bf), data = d)),
                  error = function(e) e)
  if (inherits(fit, "error")) {cat("frm ERROR:", conditionMessage(fit), "\n"); next}
  ds <- tryCatch(quiet(frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 2)),
                 error = function(e) e)
  if (inherits(ds, "error")) {cat("frm_sample ERROR:", conditionMessage(ds), "\n"); next}
  set.seed(1)
  pj <- tryCatch(quiet(posterior_predict(ds)), error = function(e) e)
  set.seed(1)
  p1 <- tryCatch(quiet(posterior_predict(ds, resp = names(fit$spec$responses)[1])),
                 error = function(e) e)
  if (inherits(pj, "error")) {cat("pp joint ERROR:", conditionMessage(pj), "\n"); next}
  if (inherits(p1, "error")) {cat("pp resp1 ERROR:", conditionMessage(p1), "\n"); next}
  cat("joint dim", dim(pj), " resp1 dim", dim(p1), "\n")
  cat("layer1: min", signif(min(pj[, , 1]), 4), " resp1 min", signif(min(p1), 4),
      "; mean", signif(mean(pj[, , 1]), 4), "vs", signif(mean(p1), 4),
      "; sd", signif(sd(pj[, , 1]), 4), "vs", signif(sd(p1), 4), "\n")
  e <- pj - quiet(posterior_epred(ds))
  cat("mean row cor of joint residuals:", round(mean(vapply(seq_len(dim(e)[2]),
      function(i) cor(e[, i, 1], e[, i, 2]), 1)), 3), "\n")
  core <- tryCatch(quiet(predict(fit, summary = FALSE, ndraws = 500)), error = function(e) e)
  if (inherits(core, "error")) cat("core predict ERROR:", conditionMessage(core), "\n") else
    cat("core predict layer1 min", signif(min(core[, , 1]), 4), "\n")
}
