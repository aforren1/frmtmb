# What brms 2.23.0 names a smooth term and its conditional_smooths()
# key, read from brms's own internals without compiling anything.
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
set.seed(1)
d <- data.frame(x = runif(50), z = runif(50),
                f = factor(sample(c("a", "b"), 50, TRUE)),
                y = rnorm(50), y2 = rnorm(50))
show <- function(form) {
  bt <- brms:::brmsterms(form)
  parts <- if (brms:::is.mvbrmsterms(bt)) {
    unlist(lapply(bt$terms, function(t) c(t$dpars, t$nlpars)),
           recursive = FALSE)
  } else {
    c(bt$dpars, bt$nlpars)
  }
  for (p in parts) {
    if (!inherits(p, "btl")) next
    sm <- tryCatch(brms:::frame_sm(p, d), error = function(e) NULL)
    if (is.null(sm) || !nrow(sm)) next
    cat(sprintf("%-40s key prefix '%s' terms: %s\n",
                deparse1(p$formula), brms:::combine_prefix(p, keep_mu = TRUE),
                paste(unique(brms:::rm_wsp(sm$term)), collapse = " | ")))
  }
}
show(bf(y ~ s(x) + s(z, by = f) + t2(x, z)))
show(bf(y ~ s(x, k = 5), sigma ~ s(z)))
show(bf(y ~ a + b, a ~ s(x), b ~ 1, nl = TRUE))
show(bf(y ~ s(x)) + bf(y2 ~ s(z)) + set_rescor(FALSE))
show(bf(y ~ s(x, f, bs = "fs")))
