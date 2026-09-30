# Reviewer re-check (lane ceplot, punch round 1): "old_levels" on a
# multivariate model whose two responses both read (1 | g), fitted() and
# predict() in one call, against brms 2.23.0 at seeds 1 to 12 (section
# (g) of dev/ceplot-rev-oldlevels3.R, whose decoder did not read an mv
# ranef()). Data seed 3.
#   Rscript dev/ceplot-rev-oldlevels4.R > dev/ceplot-rev-log/p1/oldlevels4.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(brms)
  library(frmtmb)
})
options(mc.cores = 1)
say <- function(...) cat(sprintf(...), "\n", sep = "")
set.seed(3)
ng <- 12L
d <- data.frame(x = rnorm(360), g = factor(rep(seq_len(ng), 30)))
u <- rnorm(ng, 0, 2)
d$y <- rnorm(360, 1 + 0.5 * d$x + u[d$g])
d$y2 <- rnorm(360, -1 + 0.3 * d$x + u[d$g] + rnorm(ng)[d$g])
fv <- frm(bf(y ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)), family = gaussian(),
          data = d)
cat("blocks:", vapply(fv$frame$re_blocks, `[[`, "", "term_label"), "\n")
bv <- suppressMessages(suppressWarnings(brm(
  brms::bf(y ~ x + (1 | g)) + brms::bf(y2 ~ x + (1 | g)) +
    brms::set_rescor(FALSE), data = d, family = stats::gaussian(),
  algorithm = "fixed_param", chains = 40, iter = 1, warmup = 0,
  refresh = 0, seed = 1, silent = 2)))
nd <- data.frame(x = 0, g = factor("n1"))
one <- function(k) data.frame(x = 0, g = factor(k, levels = 1:ng))
# frmtmb: per response, the seen level whose fitted value it equals
seen_f <- lapply(seq_len(ng), function(k) fitted(fv, newdata = one(k)))
lev_f <- function(f, r) {
  k <- which(vapply(seen_f, function(a) {
    identical(unname(a[1, "Estimate", r]), unname(f[1, "Estimate", r]))
  }, NA))
  if (length(k) == 1L) k else NA_integer_
}
seen_b <- lapply(seq_len(ng), function(k) {
  fitted(bv, newdata = one(k), summary = FALSE)
})
lev_b <- function(b, r) {
  k <- which(vapply(seen_b, function(a) isTRUE(all.equal(a[, 1, r], b[, 1, r],
                                                         tolerance = 0)),
                    NA))
  if (length(k) == 1L) k else NA_integer_
}
res <- t(vapply(1:12, function(s) {
  set.seed(s)
  f <- fitted(fv, newdata = nd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels")
  set.seed(s)
  b <- fitted(bv, newdata = nd, allow_new_levels = TRUE,
              sample_new_levels = "old_levels", summary = FALSE)
  c(lev_f(f, 1), lev_f(f, 2), lev_b(b, 1), lev_b(b, 2))
}, integer(4)))
say("fitted  frmtmb y/y2: %s", paste(res[, 1], res[, 2], sep = "/",
                                     collapse = " "))
say("fitted  brms   y/y2: %s", paste(res[, 3], res[, 4], sep = "/",
                                     collapse = " "))
say("one seen group for both responses: frmtmb %d of 12, brms %d of 12",
    sum(res[, 1] == res[, 2], na.rm = TRUE),
    sum(res[, 3] == res[, 4], na.rm = TRUE))
say("frmtmb's y level equals brms's at %d of 12", sum(res[, 1] == res[, 3],
                                                      na.rm = TRUE))
# predict(): the choice predict_new_level_spec() records, per response
pp <- t(vapply(1:12, function(s) {
  set.seed(s)
  a <- frmtmb:::predict_new_level_spec(fv, fv$spec$responses[["y"]], nd,
                                       NULL, TRUE, "old_levels")
  set.seed(s)
  b <- frmtmb:::predict_new_level_spec(fv, fv$spec$responses[["y2"]], nd,
                                       NULL, TRUE, "old_levels")
  c(paste(unlist(attr(a, "old_pick")), collapse = ","),
    paste(unlist(attr(b, "old_pick")), collapse = ","))
}, character(2)))
say("predict_new_level_spec y/y2 per seed (each called at the seed): %s",
    paste(pp[, 1], pp[, 2], sep = "/", collapse = " "))
set.seed(1)
pr <- tryCatch(predict(fv, newdata = nd, allow_new_levels = TRUE,
                       sample_new_levels = "old_levels", summary = FALSE,
                       ndraws = 5),
               error = function(e) conditionMessage(e))
say("predict(mv, old_levels) answers: %s", if (is.character(pr)) pr else
  paste(dim(pr), collapse = "x"))
say("done")
