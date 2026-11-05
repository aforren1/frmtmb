# Sweep, methods half: refused method rows with a brms counterpart, and
# specs corrected from the first sweep. Seed 4040.
.libPaths(c("C:/Users/adf44/source/r/wt-refuse-lib", "C:/Users/adf44/source/r/rellib-r16", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(4040); n <- 120
d <- data.frame(x = rnorm(n)); d$y <- 1 + d$x + rnorm(n); d$y2 <- d$x + rnorm(n)
d$yo <- findInterval(d$x + rnorm(n), c(-1, 0, 1)) + 1L
try1 <- function(lab, expr) { r <- tryCatch({ force(expr); "ok" }, error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 170))); cat(sprintf("%-34s %s\n", lab, r)) }
mv <- frm(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE), data = d)
rc <- frm(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = TRUE), data = d)
for (m in c("fitted", "predict", "residuals", "simulate", "conditional_effects", "emmeans")) {
  fn <- switch(m, emmeans = function(f, ...) emmeans::emmeans(f, ~ x, ...), get(m))
  try1(paste0("mvbf x ", m), fn(mv))
  try1(paste0("rescor x ", m), fn(rc))
}
try1("thres(3) x mvbf", frm(mvbf(bf(yo | thres(3) ~ x, family = cumulative()), bf(y ~ x), rescor = FALSE), data = d))
try1("sum_to_zero x mvbf", frm(mvbf(bf(yo ~ x, family = cumulative(threshold = "sum_to_zero")), bf(y ~ x), rescor = FALSE), data = d))
b <- tryCatch(brms::stancode(brms::bf(yo | thres(3) ~ x, family = brms::cumulative()) + brms::bf(y ~ x) + brms::set_rescor(FALSE), data = d), error = function(e) e)
cat("brms thres(3) x mvbf:", if (inherits(b, "error")) conditionMessage(b) else "builds", "\n")
