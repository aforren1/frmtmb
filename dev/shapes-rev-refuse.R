# Reviewer, priority 7: every refusal the lane left in place. For each
# one: does it fire, is the condition CLASSED, does the message name a
# reason, and is the same question reachable through another spelling
# that ANSWERS instead?
#
#   Rscript dev/shapes-rev-refuse.R

.libPaths(c("C:/Users/adf44/source/r/shapes-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
suppressPackageStartupMessages(library(frmtmb.sample))
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
source(file.path(TREE, "dev/shapes-rev-fixtures.R"))
dd <- rev_data()
dd$y2 <- with(dd, rnorm(nrow(dd), 0.3 - 0.2 * x, 1))

g  <- frm(bf(y ~ x + f) + gaussian(), data = dd)
mx <- frm(bf(ymix ~ x + (1 | g)) + gaussian(), data = dd)
mv <- frm(bf(mvbind(y, y2) ~ x) + gaussian(), data = dd)

probe <- function(lab, expr) {
  r <- tryCatch(list(v = suppressWarnings(expr), c = NULL),
                error = function(c) list(v = NULL, c = c))
  if (is.null(r$c)) {
    cat(sprintf("ANSWERS  %-46s -> %s\n", lab,
                paste(class(r$v)[1], paste(dim(r$v) %||% length(r$v),
                                           collapse = "x"))))
    return(invisible("answers"))
  }
  cl <- setdiff(class(r$c), c("condition"))
  cat(sprintf("REFUSES  %-46s classes: %s\n", lab,
              paste(cl, collapse = "/")))
  cat("         ", substr(conditionMessage(r$c), 1, 220), "\n")
  invisible("refuses")
}
`%||%` <- function(a, b) if (is.null(a)) b else a

cat("======== the refusals the lane LEFT, by name ========\n")
probe("fitted(ndraws = 10)", fitted(g, ndraws = 10))
probe("fitted(draw_ids = 1:3)", fitted(g, draw_ids = 1:3))
probe("fitted(summary = FALSE)", fitted(g, summary = FALSE))
probe("residuals(ndraws = 10)", residuals(g, ndraws = 10))
probe("residuals(newdata = dd)", residuals(g, newdata = dd))
probe("fitted(mv)  [multivariate]", fitted(mv))
probe("residuals(mv)  [multivariate]", residuals(mv))
probe("predict(mx, sample_new_levels = 'uncertainty')",
      predict(mx, sample_new_levels = "uncertainty"))
probe("predict(mx, sample_new_levels = FALSE)",
      predict(mx, sample_new_levels = FALSE))
probe("predict(g, draw_ids = 1:3)", predict(g, draw_ids = 1:3))
probe("predict(g, sort = TRUE)", predict(g, sort = TRUE))
probe("predict(g, cores = 2)", predict(g, cores = 2))
probe("predict(g, type = 'response')", predict(g, type = "response"))
probe("predict(g, se.fit = TRUE)", predict(g, se.fit = TRUE))
probe("predict(g, dpar = 'mu')", predict(g, dpar = "mu"))
probe("predict(g, scale = 'response')", predict(g, scale = "response"))
probe("predict(g, re.form = NA)", predict(g, re.form = NA))
probe("predict(g, allow.new.levels = TRUE)",
      predict(g, allow.new.levels = TRUE))
probe("predict(mv)  [multivariate]", predict(mv))

cat("\n======== is the refused question reachable another way? ========\n")
# multivariate fitted/residuals: does another spelling answer?
probe("frm_linpred(mv, type='response')",
      frm_linpred(mv, type = "response"))
probe("frm_linpred(mv, resp='y2', type='response')",
      frm_linpred(mv, resp = "y2", type = "response"))
probe("predict(mv, resp = 'y2')", predict(mv, resp = "y2"))
probe("fitted(mv, resp = 'y2')", fitted(mv, resp = "y2"))
probe("residuals(mv, resp = 'y2')", residuals(mv, resp = "y2"))
probe("conditional_effects(mv)", conditional_effects(mv))
# residuals(newdata): reachable through fitted(newdata)?
probe("fitted(g, newdata = head(dd, 5))", fitted(g, newdata = head(dd, 5)))
# fitted(ndraws): reachable through predict(ndraws)?
probe("predict(g, ndraws = 10)", predict(g, ndraws = 10))
# sample_new_levels: reachable on draws?
cat("\n-- frmtmb.sample draws --\n")
dr <- tryCatch(suppressWarnings(frm_sample(mx, chains = 1, iter = 300,
                                           refresh = 0, seed = 1)),
               error = function(e) { cat("  frm_sample failed: ",
                                         conditionMessage(e), "\n"); NULL })
if (!is.null(dr)) {
  probe("predict(draws, sample_new_levels = 'uncertainty')",
        predict(dr, sample_new_levels = "uncertainty"))
  probe("posterior_predict(draws, sample_new_levels = 'uncertainty')",
        posterior_predict(dr, sample_new_levels = "uncertainty"))
  probe("fitted(draws, ndraws = 10)", fitted(dr, ndraws = 10))
  probe("residuals(draws, newdata = dd)", residuals(dr, newdata = dd))
  probe("fitted(draws)", fitted(dr))
  probe("predict(draws)", predict(dr))
  probe("residuals(draws)", residuals(dr))
}
