# What brms 2.23.0 does with pp_check(newdata =) and with re_formula in
# posterior_predict(), on its own example fit (no Stan compile needed:
# posterior_predict() runs in R on the stored draws).
#   Rscript dev/simnewdata-brms.R > dev/simnewdata-log/brms.txt 2>&1
source("dev/simnewdata-prelude.R")
suppressMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
fit1 <- brms:::rename_pars(brms:::brmsfit_example1)
cat("formula:", deparse1(formula(fit1)$formula), "\n")
nd <- fit1$data[1:10, ]
try_ <- function(expr) {
  tryCatch(expr, error = function(e) {
    cat("  ERROR [", paste(class(e), collapse = "/"), "]: ",
        conditionMessage(e), "\n", sep = "")
    invisible(NULL)
  })
}

cat("\n== :675 pp_check(fit1, newdata = fit1$data[1:10, ]) ==\n")
set.seed(1)
p <- try_(pp_check(fit1, newdata = nd, ndraws = 10))
cat("  ggplot:", inherits(p, "ggplot"), "\n")
cat("  observed series length:", sum(p$data$is_y), "\n")

cat("\n== :680/:682 violin_grouped, group read off newdata ==\n")
set.seed(1)
p2 <- try_(pp_check(fit1, type = "violin_grouped", group = "visit",
                    newdata = nd, ndraws = 10))
cat("  ggplot:", inherits(p2, "ggplot"), "\n")
cat("  groups in plot data:", paste(sort(unique(as.character(p2$data$group))),
                                  collapse = ","), "\n")
cat("  newdata visit levels present:",
    paste(sort(unique(as.character(nd$visit))), collapse = ","), "\n")

cat("\n== a group column absent from newdata ==\n")
nd2 <- nd
nd2$visit <- NULL
try_(pp_check(fit1, type = "violin_grouped", group = "visit",
              newdata = nd2, ndraws = 10, re_formula = NA))

cat("\n== the response absent from newdata ==\n")
nd3 <- nd
nd3$count <- NULL
r <- try_(pp_check(fit1, newdata = nd3, ndraws = 10))
cat("  returned:", !is.null(r), "\n")
r <- try_(pp_check(fit1, newdata = nd3, ndraws = 10, prefix = "ppd"))
cat("  ppd returned:", !is.null(r), "\n")

cat("\n== an unseen grouping level in newdata ==\n")
nd4 <- nd
nd4$visit <- factor(rep("99", 10))
try_(posterior_predict(fit1, newdata = nd4, ndraws = 5))
r <- try_(posterior_predict(fit1, newdata = nd4, ndraws = 5,
                            allow_new_levels = TRUE))
cat("  allow_new_levels = TRUE returned:", !is.null(r), "\n")

cat("\n== re_formula in posterior_predict: ~1 against NA ==\n")
pp_at <- function(rf) {
  set.seed(7)
  posterior_predict(fit1, newdata = nd, re_formula = rf, draw_ids = seq_len(ndraws(fit1)))
}
cat("  ~1 identical to NA:", identical(pp_at(~1), pp_at(NA)), "\n")
cat("  ~0 identical to NA:", identical(pp_at(~0), pp_at(NA)), "\n")
cat("  NULL identical to NA:", identical(pp_at(NULL), pp_at(NA)), "\n")


cat("\n== a partial re_formula keeps the named columns ==\n")
lp_full <- posterior_linpred(fit1, newdata = nd, draw_ids = seq_len(ndraws(fit1)))
lp_int <- posterior_linpred(fit1, newdata = nd, draw_ids = seq_len(ndraws(fit1)),
                            re_formula = ~ (1 | visit))
lp_na <- posterior_linpred(fit1, newdata = nd, draw_ids = seq_len(ndraws(fit1)),
                           re_formula = NA)
cat("  ~(1 | visit) differs from NULL:", !isTRUE(all.equal(lp_int, lp_full)),
    "and from NA:", !isTRUE(all.equal(lp_int, lp_na)), "\n")
r <- try_(posterior_linpred(fit1, newdata = nd, draw_ids = seq_len(ndraws(fit1)),
                            re_formula = ~ (1 | nosuch)))
cat("  ~(1 | nosuch) equals NA:", isTRUE(all.equal(r, lp_na)), "\n")

cat("\n== multinomial: pp_check refuses outright ==\n")
set.seed(3)
N <- 30
dm <- data.frame(x = rnorm(N), n = 10L)
dm$Y <- t(sapply(seq_len(N), function(i) rmultinom(1, 10, c(.2, .3, .5))))
colnames(dm$Y) <- c("a", "b", "c")
fm <- try_(brm(bf(Y | trials(n) ~ x), data = dm, family = multinomial(),
               chains = 0, silent = 2))
if (!is.null(fm)) {
  try_(pp_check(fm, type = "error_binned"))
  try_(pp_check(fm))
}
